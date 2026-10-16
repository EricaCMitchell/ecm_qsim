#include "ecm_qsim/state_vector.hpp"

#include <unistd.h>

#include <algorithm>
#include <bit>
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
