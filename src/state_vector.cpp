#include "ecm_qsim/state_vector.hpp"

#include <unistd.h>

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <iomanip>
#include <limits>
#include <new>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>

#include "ecm_qsim/circuit.hpp"

namespace ecm_qsim {

namespace {

static_assert(sizeof(std::complex<double>) == 16);

// Probabilities below this are rounding noise; to_string() omits them.
constexpr double kNegligibleProbability = 1e-12;

// Bytes per amplitude that a simulation holds while it samples: the
// amplitude (16) plus its entry in one probabilities() copy (8).
constexpr std::uint64_t kWorkingBytesPerAmplitude = 24;

void check_qubit_count(int n_qubits, const std::string& where) {
  if (n_qubits < 1 || n_qubits > kMaxQubits) {
    throw std::invalid_argument(where + ": qubit count must be 1.." +
                                std::to_string(kMaxQubits) + ", got " +
                                std::to_string(n_qubits));
  }
}

// Physical memory in bytes, or the largest uint64_t if the OS does not
// report it (which effectively turns off the size check).
std::uint64_t physical_memory_bytes() {
  const long pages = sysconf(_SC_PHYS_PAGES);
  const long page_size = sysconf(_SC_PAGESIZE);
  if (pages <= 0 || page_size <= 0) {
    return std::numeric_limits<std::uint64_t>::max();
  }
  return static_cast<std::uint64_t>(pages) *
         static_cast<std::uint64_t>(page_size);
}

// |0...0> on n_qubits qubits. Checks the size before allocating, so callers
// receive std::length_error instead of std::bad_alloc or, on an OS that
// overcommits memory, a kill signal.
Eigen::VectorXcd zero_state(int n_qubits) {
  check_qubit_count(n_qubits, "StateVector");
  const int limit = StateVector::max_qubits();
  if (n_qubits > limit) {
    throw std::length_error(
        "StateVector: " + std::to_string(n_qubits) + " qubits need 2^" +
        std::to_string(n_qubits) + " amplitudes; at " +
        std::to_string(kWorkingBytesPerAmplitude) + " bytes each, at most " +
        std::to_string(limit) + " qubits fit in half of physical memory");
  }
  // Backstop only: macOS overcommits, so an oversized allocation usually
  // succeeds and the OS kills the process later instead of throwing.
  try {
    Eigen::VectorXcd amplitudes =
        Eigen::VectorXcd::Zero(Eigen::Index{1} << n_qubits);
    amplitudes(0) = 1.0;
    return amplitudes;
  } catch (const std::bad_alloc&) {
    throw std::length_error("StateVector: could not allocate 2^" +
                            std::to_string(n_qubits) + " amplitudes for " +
                            std::to_string(n_qubits) + " qubits");
  }
}

int qubits_for_size(Eigen::Index size) {
  if (size < 2 || !std::has_single_bit(static_cast<std::uint64_t>(size))) {
    throw std::invalid_argument(
        "StateVector: amplitude count must be 2^n with n >= 1, got " +
        std::to_string(size));
  }
  return std::countr_zero(static_cast<std::uint64_t>(size));
}

void check_qubit(int qubit, int n_qubits, const char* where) {
  if (qubit < 0 || qubit >= n_qubits) {
    throw std::out_of_range(std::string(where) + ": qubit " +
                            std::to_string(qubit) + " out of range for " +
                            std::to_string(n_qubits) + " qubits");
  }
}

}  // namespace

StateVector::StateVector(int n_qubits)
    : n_qubits_(n_qubits), amplitudes_(zero_state(n_qubits)) {}

StateVector::StateVector(Eigen::VectorXcd amplitudes)
    : n_qubits_(qubits_for_size(amplitudes.size())),
      amplitudes_(std::move(amplitudes)) {}

int StateVector::max_qubits() {
  // Divide the budget instead of multiplying 2^n by 24, so that n up to 63
  // cannot overflow. bit_width(x) - 1 is floor(log2(x)).
  const std::uint64_t max_amplitudes =
      physical_memory_bytes() / 2 / kWorkingBytesPerAmplitude;
  const int n = static_cast<int>(std::bit_width(max_amplitudes)) - 1;
  return std::min(n, kMaxQubits);
}

std::complex<double> StateVector::amplitude(std::uint64_t index) const {
  if (index >= size()) {
    throw std::out_of_range("StateVector::amplitude: index " +
                            std::to_string(index) + " out of range for " +
                            std::to_string(n_qubits_) + " qubits");
  }
  return amplitudes_(static_cast<Eigen::Index>(index));
}

Eigen::VectorXcd StateVector::to_vector() const {
  return amplitudes_;
}

Eigen::VectorXd StateVector::probabilities() const {
  return amplitudes_.cwiseAbs2();
}

void StateVector::apply_1q(const Eigen::Matrix2cd& u, int k) {
  check_qubit(k, n_qubits_, "apply_1q");
  const std::uint64_t bit = std::uint64_t{1} << k;
  const std::uint64_t dim = size();
  for (std::uint64_t i = 0; i < dim; ++i) {
    if ((i & bit) != 0) {
      continue;
    }
    const auto i0 = static_cast<Eigen::Index>(i);
    const auto i1 = static_cast<Eigen::Index>(i + bit);
    const Eigen::Vector2cd out =
        u * Eigen::Vector2cd(amplitudes_(i0), amplitudes_(i1));
    amplitudes_(i0) = out(0);
    amplitudes_(i1) = out(1);
  }
}

void StateVector::apply_2q(const Eigen::Matrix4cd& m, int j, int k) {
  check_qubit(j, n_qubits_, "apply_2q");
  check_qubit(k, n_qubits_, "apply_2q");
  if (j == k) {
    throw std::invalid_argument("apply_2q: both qubits are " +
                                std::to_string(j));
  }
  const std::uint64_t bit_j = std::uint64_t{1} << j;
  const std::uint64_t bit_k = std::uint64_t{1} << k;
  const std::uint64_t dim = size();
  for (std::uint64_t i = 0; i < dim; ++i) {
    if ((i & (bit_j | bit_k)) != 0) {
      continue;
    }
    // m's row order: (qubit j, qubit k) = 00, 01, 10, 11.
    const std::array<Eigen::Index, 4> index = {
        static_cast<Eigen::Index>(i),
        static_cast<Eigen::Index>(i + bit_k),
        static_cast<Eigen::Index>(i + bit_j),
        static_cast<Eigen::Index>(i + bit_j + bit_k),
    };
    const Eigen::Vector4cd out =
        m * Eigen::Vector4cd(amplitudes_(index[0]), amplitudes_(index[1]),
                             amplitudes_(index[2]), amplitudes_(index[3]));
    for (int r = 0; r < 4; ++r) {
      amplitudes_(index[r]) = out(r);
    }
  }
}

void StateVector::apply_ry(double theta, int k) {
  check_qubit(k, n_qubits_, "apply_ry");
  const double c = std::cos(theta / 2);
  const double s = std::sin(theta / 2);
  const std::uint64_t bit = std::uint64_t{1} << k;
  const std::uint64_t dim = size();
  // Blocks of 2 * bit indices: the first half has bit k clear, and its
  // partners are the second half.
  for (std::uint64_t block = 0; block < dim; block += 2 * bit) {
    for (std::uint64_t i = block; i < block + bit; ++i) {
      const auto i0 = static_cast<Eigen::Index>(i);
      const auto i1 = static_cast<Eigen::Index>(i + bit);
      const std::complex<double> a0 = amplitudes_(i0);
      const std::complex<double> a1 = amplitudes_(i1);
      amplitudes_(i0) = c * a0 - s * a1;
      amplitudes_(i1) = s * a0 + c * a1;
    }
  }
}

void StateVector::apply_cnot(int control, int target) {
  check_qubit(control, n_qubits_, "apply_cnot");
  check_qubit(target, n_qubits_, "apply_cnot");
  if (control == target) {
    throw std::invalid_argument("apply_cnot: both qubits are " +
                                std::to_string(control));
  }
  const std::uint64_t control_bit = std::uint64_t{1} << control;
  const std::uint64_t target_bit = std::uint64_t{1} << target;
  const std::uint64_t low_bit = std::min(control_bit, target_bit);
  const std::uint64_t high_bit = std::max(control_bit, target_bit);
  const std::uint64_t dim = size();
  // high_base + middle + offset runs over every index with both gate bits
  // clear, each once: 2^(n-2) indices.
  for (std::uint64_t high_base = 0; high_base < dim; high_base += 2 * high_bit) {
    for (std::uint64_t middle = 0; middle < high_bit; middle += 2 * low_bit) {
      for (std::uint64_t offset = 0; offset < low_bit; ++offset) {
        const std::uint64_t i = high_base + middle + offset + control_bit;
        std::swap(amplitudes_(static_cast<Eigen::Index>(i)),
                  amplitudes_(static_cast<Eigen::Index>(i + target_bit)));
      }
    }
  }
}

std::string StateVector::to_string() const {
  std::ostringstream out;
  out << std::fixed << std::setprecision(4);
  for (Eigen::Index i = 0; i < amplitudes_.size(); ++i) {
    const double p = std::norm(amplitudes_(i));
    if (p < kNegligibleProbability) {
      continue;
    }
    out << '|' << to_bitstring(static_cast<std::uint64_t>(i), n_qubits_)
        << ">  p=" << p << '\n';
  }
  return out.str();
}

std::ostream& operator<<(std::ostream& os, const StateVector& state) {
  return os << state.to_string();
}

std::string to_bitstring(std::uint64_t index, int n_qubits) {
  check_qubit_count(n_qubits, "to_bitstring");
  if ((index >> n_qubits) != 0) {
    throw std::out_of_range("to_bitstring: index " + std::to_string(index) +
                            " needs more than " + std::to_string(n_qubits) +
                            " bits");
  }
  std::string bits(static_cast<std::size_t>(n_qubits), '0');
  for (int q = 0; q < n_qubits; ++q) {
    if (((index >> q) & 1) != 0) {
      bits[static_cast<std::size_t>(n_qubits - 1 - q)] = '1';
    }
  }
  return bits;
}

}  // namespace ecm_qsim
