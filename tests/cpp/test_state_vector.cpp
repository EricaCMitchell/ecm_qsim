#include <unistd.h>

#include <Eigen/Core>
#include <array>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_exception.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <initializer_list>
#include <numbers>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unsupported/Eigen/KroneckerProduct>
#include <vector>

#include "ecm_qsim/circuit.hpp"
#include "ecm_qsim/gates.hpp"
#include "ecm_qsim/state_vector.hpp"
#include "eigen_helpers.hpp"

using namespace ecm_qsim;
using Catch::Matchers::ContainsSubstring;

namespace {

using cd = std::complex<double>;

Eigen::VectorXcd vec(std::initializer_list<cd> entries) {
  Eigen::VectorXcd v(static_cast<Eigen::Index>(entries.size()));
  Eigen::Index i = 0;
  for (const cd& x : entries) {
    v(i++) = x;
  }
  return v;
}

const Eigen::Matrix2cd kId2 = Eigen::Matrix2cd::Identity();

// The 2^n x 2^n operator that acts as factors[q] on qubit q.
// kroneckerProduct(A, B) puts A on the high bit, so the product runs from
// qubit n-1 down to qubit 0. Test-only oracle: O(4^n) memory.
Eigen::MatrixXcd kron_all(const std::vector<Eigen::Matrix2cd>& factors) {
  Eigen::MatrixXcd result = Eigen::MatrixXcd::Identity(1, 1);
  for (auto it = factors.rbegin(); it != factors.rend(); ++it) {
    result = Eigen::kroneckerProduct(result, *it).eval();
  }
  return result;
}

// kron(I, ..., u, ..., I) with u on qubit k of n.
Eigen::MatrixXcd embed_1q(const Eigen::Matrix2cd& u, int k, int n) {
  std::vector<Eigen::Matrix2cd> factors(static_cast<std::size_t>(n), kId2);
  factors[static_cast<std::size_t>(k)] = u;
  return kron_all(factors);
}

// Every test RNG starts from this seed, so failures reproduce.
constexpr std::uint64_t kSeed = 42;

// Entries with independent standard-normal real and imaginary parts. The
// kernels are linear, so a general matrix tests them more strictly than a
// gate. Gate unitarity is tested in test_gates.cpp.
Eigen::MatrixXcd random_matrix(Eigen::Index rows, Eigen::Index cols,
                               std::mt19937_64& rng) {
  std::normal_distribution<double> normal;
  Eigen::MatrixXcd m(rows, cols);
  for (Eigen::Index r = 0; r < rows; ++r) {
    for (Eigen::Index c = 0; c < cols; ++c) {
      m(r, c) = cd(normal(rng), normal(rng));
    }
  }
  return m;
}

Eigen::VectorXcd random_state(int n, std::mt19937_64& rng) {
  return random_matrix(Eigen::Index{1} << n, 1, rng);
}

// |a><b| on one qubit.
Eigen::Matrix2cd outer(int a, int b) {
  Eigen::Matrix2cd e = Eigen::Matrix2cd::Zero();
  e(a, b) = 1.0;
  return e;
}

// The 2^n x 2^n operator of m on qubits (j, k), built without the kernel's
// index arithmetic: the sum over entries m(r, c) of |r_j><c_j| on qubit j
// (x) |r_k><c_k| on qubit k, identity elsewhere. r_j is the high bit of r,
// matching matrix_2q's row order (first-listed qubit, second-listed qubit).
Eigen::MatrixXcd embed_2q(const Eigen::Matrix4cd& m, int j, int k, int n) {
  const Eigen::Index dim = Eigen::Index{1} << n;
  Eigen::MatrixXcd result = Eigen::MatrixXcd::Zero(dim, dim);
  for (int r = 0; r < 4; ++r) {
    for (int c = 0; c < 4; ++c) {
      std::vector<Eigen::Matrix2cd> factors(static_cast<std::size_t>(n), kId2);
      factors[static_cast<std::size_t>(j)] = outer(r >> 1, c >> 1);
      factors[static_cast<std::size_t>(k)] = outer(r & 1, c & 1);
      result += m(r, c) * kron_all(factors);
    }
  }
  return result;
}

}  // namespace

TEST_CASE("A new StateVector is the all-zero basis state") {
  const StateVector state(3);
  REQUIRE(state.n_qubits() == 3);
  REQUIRE(state.size() == 8);
  require_close(state.to_vector(), Eigen::VectorXcd::Unit(8, 0));
}

TEST_CASE("StateVector rejects qubit counts outside 1..kMaxQubits") {
  const int n = GENERATE(-1, 0, 64);
  CAPTURE(n);
  REQUIRE_THROWS_AS(StateVector(n), std::invalid_argument);
}

TEST_CASE("StateVector refuses sizes beyond its memory limit") {
  // 40 qubits need 2^40 * 16 bytes (16 TiB). For 63 qubits the byte count
  // does not fit in 64 bits.
  const int n = GENERATE(40, 63);
  CAPTURE(n);
  REQUIRE_THROWS_AS(StateVector(n), std::length_error);
  REQUIRE_THROWS_WITH(StateVector(n),
                      ContainsSubstring(std::to_string(n) + " qubits need 2^" +
                                        std::to_string(n) + " amplitudes"));
}

TEST_CASE("max_qubits fits 16 bytes per amplitude in memory less 1 GiB") {
  // 2^limit amplitudes fit in the budget and 2^(limit + 1) do not. The
  // state at the limit is not built: it can take almost all of physical
  // memory.
  const int limit = StateVector::max_qubits();
  CAPTURE(limit);
  REQUIRE(limit >= 20);  // 2^20 * 16 bytes is 16 MiB
  REQUIRE(limit < kMaxQubits);
  const std::uint64_t budget =
      static_cast<std::uint64_t>(sysconf(_SC_PHYS_PAGES)) *
          static_cast<std::uint64_t>(sysconf(_SC_PAGESIZE)) -
      (std::uint64_t{1} << 30);
  REQUIRE((std::uint64_t{16} << limit) <= budget);
  REQUIRE((std::uint64_t{16} << (limit + 1)) > budget);
  // The constructor refuses limit + 1 before it allocates.
  REQUIRE_THROWS_MATCHES(
      StateVector(limit + 1), std::length_error,
      Catch::Matchers::MessageMatches(Catch::Matchers::ContainsSubstring(
          "at 16 bytes each, at most " + std::to_string(limit) +
          " qubits fit in physical memory less 1 GiB")));
}

TEST_CASE("StateVector takes amplitudes as given") {
  const Eigen::VectorXcd amplitudes = vec({0.6, 0, 0, cd(0, 0.8)});
  const StateVector state(amplitudes);
  REQUIRE(state.n_qubits() == 2);
  REQUIRE(state.size() == 4);
  REQUIRE(state.amplitude(3) == cd(0, 0.8));
  require_close(state.to_vector(), amplitudes);
}

TEST_CASE("StateVector rejects amplitude counts that are not 2^n, n >= 1") {
  const int size = GENERATE(0, 1, 3, 6);
  CAPTURE(size);
  const Eigen::VectorXcd zeros = Eigen::VectorXcd::Zero(size);
  REQUIRE_THROWS_AS(StateVector(zeros), std::invalid_argument);
}

TEST_CASE("amplitude rejects an index past the end") {
  const StateVector state(2);
  REQUIRE_THROWS_AS(state.amplitude(4), std::out_of_range);
}

TEST_CASE("to_vector returns a copy") {
  const StateVector state(1);
  Eigen::VectorXcd copy = state.to_vector();
  copy(0) = 0.0;
  REQUIRE(state.amplitude(0) == cd(1, 0));
}

TEST_CASE("probabilities are squared magnitudes") {
  const StateVector state(vec({0.6, 0, 0, cd(0, -0.8)}));
  const Eigen::VectorXd p = state.probabilities();
  REQUIRE(p.size() == 4);
  REQUIRE(p(0) == near(0.36));
  REQUIRE(p(1) == near(0.0));
  REQUIRE(p(2) == near(0.0));
  REQUIRE(p(3) == near(0.64));
}

TEST_CASE("to_bitstring puts qubit 0 rightmost") {
  REQUIRE(to_bitstring(1, 3) == "001");
  REQUIRE(to_bitstring(3, 3) == "011");
  REQUIRE(to_bitstring(4, 3) == "100");
  REQUIRE(to_bitstring(5, 3) == "101");
  REQUIRE(to_bitstring(6, 3) == "110");
  REQUIRE(to_bitstring(1, 1) == "1");
}

TEST_CASE("to_bitstring handles 63 qubits") {
  REQUIRE(to_bitstring(std::uint64_t{1} << 62, 63) ==
          "1" + std::string(62, '0'));
  REQUIRE(to_bitstring((std::uint64_t{1} << 63) - 1, 63) ==
          std::string(63, '1'));
}

TEST_CASE("to_bitstring rejects bad arguments") {
  REQUIRE_THROWS_AS(to_bitstring(8, 3), std::out_of_range);
  REQUIRE_THROWS_AS(to_bitstring(0, 0), std::invalid_argument);
  REQUIRE_THROWS_AS(to_bitstring(0, 64), std::invalid_argument);
}

TEST_CASE("to_string lists non-negligible basis states in index order") {
  const double r = 1 / std::sqrt(2.0);
  // Index 1 has p = 1e-18: rounding noise, so it is omitted.
  const StateVector state(vec({r, 1e-9, 0, 0, 0, 0, 0, -r}));
  REQUIRE(state.to_string() == "|000>  p=0.5000\n|111>  p=0.5000\n");
}

TEST_CASE("to_string labels each state with its bitstring") {
  // Indices 1 and 5 are 001 and 101.
  const StateVector state(vec({0, 0.6, 0, 0, 0, cd(0, 0.8), 0, 0}));
  REQUIRE(state.to_string() == "|001>  p=0.3600\n|101>  p=0.6400\n");
}

TEST_CASE("to_string keeps a small but real probability") {
  // p = 1e-6 is above the 1e-12 threshold, so it is listed, although it
  // rounds to 0.0000 at 4 decimals.
  const StateVector state(vec({1, 1e-3}));
  REQUIRE(state.to_string() == "|0>  p=1.0000\n|1>  p=0.0000\n");
}

TEST_CASE("operator<< writes to_string()") {
  const StateVector state(2);
  std::ostringstream out;
  out << state;
  REQUIRE(out.str() == "|00>  p=1.0000\n");
}

TEST_CASE("Oracle: X on qubit 0 of 3 swaps neighbouring indices") {
  // Eq. 17: I (x) I (x) X swaps 0<->1, 2<->3, 4<->5, 6<->7.
  Eigen::MatrixXcd expected = Eigen::MatrixXcd::Zero(8, 8);
  for (int i = 0; i < 8; ++i) {
    expected(i, i ^ 1) = 1.0;
  }
  require_close(embed_1q(matrix_1q(GateType::X), 0, 3), expected);
}

TEST_CASE("apply_1q matches the Kronecker-product oracle") {
  std::mt19937_64 rng(kSeed);
  for (int n = 1; n <= 5; ++n) {
    for (int k = 0; k < n; ++k) {
      CAPTURE(n, k);
      const Eigen::VectorXcd psi = random_state(n, rng);
      const Eigen::Matrix2cd u = random_matrix(2, 2, rng);
      StateVector state(psi);
      state.apply_1q(u, k);
      require_close(state.to_vector(), embed_1q(u, k, n) * psi);
    }
  }
}

TEST_CASE("X on qubit 0 of 3 qubits gives bitstring 001") {
  StateVector state(3);
  state.apply_1q(matrix_1q(GateType::X), 0);
  require_close(state.to_vector(), Eigen::VectorXcd::Unit(8, 1));
  REQUIRE(state.to_string() == "|001>  p=1.0000\n");
}

TEST_CASE("X on the highest qubit sets the leftmost bit") {
  StateVector state(3);
  state.apply_1q(matrix_1q(GateType::X), 2);
  REQUIRE(state.to_string() == "|100>  p=1.0000\n");
}

TEST_CASE("to_string omits the rounding residue left by RX(pi)") {
  // cos(pi / 2) is 6e-17, not 0, so |0> keeps p of about 4e-33.
  StateVector state(1);
  state.apply_1q(matrix_1q(GateType::RX, std::numbers::pi), 0);
  REQUIRE(std::norm(state.amplitude(0)) > 0.0);
  REQUIRE(state.to_string() == "|1>  p=1.0000\n");
}

TEST_CASE("apply_1q then the adjoint restores the state") {
  const GateType gate =
      GENERATE(GateType::X, GateType::Y, GateType::Z, GateType::H, GateType::S,
               GateType::T, GateType::RX, GateType::RY, GateType::RZ);
  CAPTURE(gate);
  std::mt19937_64 rng(kSeed);
  const Eigen::VectorXcd psi = random_state(3, rng);
  const Eigen::Matrix2cd u = matrix_1q(gate, 0.7);
  for (int k = 0; k < 3; ++k) {
    CAPTURE(k);
    StateVector state(psi);
    state.apply_1q(u, k);
    state.apply_1q(u.adjoint(), k);
    require_close(state.to_vector(), psi);
  }
}

TEST_CASE("apply_1q rejects a qubit outside the state") {
  const int k = GENERATE(-1, 3, 64);
  CAPTURE(k);
  StateVector state(3);
  const Eigen::Matrix2cd x = matrix_1q(GateType::X);
  REQUIRE_THROWS_AS(state.apply_1q(x, k), std::out_of_range);
}

TEST_CASE("Oracle: on 2 qubits, embed_2q(m, 1, 0) is m itself") {
  // Qubit 1 is the high bit of the 4-entry index, as j is in m.
  std::mt19937_64 rng(kSeed);
  const Eigen::Matrix4cd m = random_matrix(4, 4, rng);
  require_close(embed_2q(m, 1, 0, 2), m);
}

TEST_CASE("Oracle: CNOT 0 2 on 3 qubits swaps 1<->5 and 3<->7") {
  Eigen::MatrixXcd expected = Eigen::MatrixXcd::Identity(8, 8);
  expected.row(1).swap(expected.row(5));
  expected.row(3).swap(expected.row(7));
  require_close(embed_2q(matrix_2q(GateType::CNOT), 0, 2, 3), expected);
}

TEST_CASE("apply_2q matches the oracle for every ordered qubit pair") {
  std::mt19937_64 rng(kSeed);
  for (int n = 2; n <= 5; ++n) {
    for (int j = 0; j < n; ++j) {
      for (int k = 0; k < n; ++k) {
        if (j == k) {
          continue;
        }
        CAPTURE(n, j, k);
        const Eigen::VectorXcd psi = random_state(n, rng);
        const Eigen::Matrix4cd m = random_matrix(4, 4, rng);
        StateVector state(psi);
        state.apply_2q(m, j, k);
        require_close(state.to_vector(), embed_2q(m, j, k, n) * psi);
      }
    }
  }
}

TEST_CASE("Worked example: CNOT 0 2 swaps psi[1]<->psi[5], psi[3]<->psi[7]") {
  StateVector state(vec({1, 2, 3, 4, 5, 6, 7, 8}));
  state.apply_2q(matrix_2q(GateType::CNOT), 0, 2);
  require_close(state.to_vector(), vec({1, 6, 3, 8, 5, 2, 7, 4}));
}

TEST_CASE("CNOT with the control above the target") {
  const Eigen::Matrix4cd cnot = matrix_2q(GateType::CNOT);
  // CNOT 2 0 on |100>: control (qubit 2) is 1, so qubit 0 flips.
  StateVector flipped(Eigen::VectorXcd::Unit(8, 4));
  flipped.apply_2q(cnot, 2, 0);
  REQUIRE(flipped.to_string() == "|101>  p=1.0000\n");
  // CNOT 2 0 on |001>: control is 0, so nothing changes.
  StateVector unchanged(Eigen::VectorXcd::Unit(8, 1));
  unchanged.apply_2q(cnot, 2, 0);
  REQUIRE(unchanged.to_string() == "|001>  p=1.0000\n");
}

TEST_CASE("H and two CNOTs make the 3-qubit GHZ state") {
  StateVector state(3);
  state.apply_1q(matrix_1q(GateType::H), 0);
  state.apply_2q(matrix_2q(GateType::CNOT), 0, 1);
  state.apply_2q(matrix_2q(GateType::CNOT), 1, 2);
  REQUIRE(state.to_string() == "|000>  p=0.5000\n|111>  p=0.5000\n");
}

TEST_CASE("apply_2q then the adjoint restores the state") {
  const GateType gate = GENERATE(GateType::CNOT, GateType::CZ, GateType::CY);
  CAPTURE(gate);
  std::mt19937_64 rng(kSeed);
  const Eigen::VectorXcd psi = random_state(3, rng);
  const Eigen::Matrix4cd m = matrix_2q(gate);
  for (int j = 0; j < 3; ++j) {
    for (int k = 0; k < 3; ++k) {
      if (j == k) {
        continue;
      }
      CAPTURE(j, k);
      StateVector state(psi);
      state.apply_2q(m, j, k);
      state.apply_2q(m.adjoint(), j, k);
      require_close(state.to_vector(), psi);
    }
  }
}

TEST_CASE("Every gate keeps the total probability within 1e-12 of 1") {
  // 500 random gates of every type on random qubits, with controls both
  // above and below their targets. probabilities() holds the values that
  // sample_index reads, so its sum is checked after each gate: a kernel that
  // loses norm fails at the gate that caused it, and rounding drift fails as
  // soon as it passes 1e-12.
  constexpr int kQubits = 6;
  constexpr int kGates = 500;
  constexpr std::array kOneQubit = {GateType::X,  GateType::Y,  GateType::Z,
                                    GateType::H,  GateType::S,  GateType::T,
                                    GateType::RX, GateType::RY, GateType::RZ};
  constexpr std::array kTwoQubit = {GateType::CNOT, GateType::CZ, GateType::CY};
  std::mt19937_64 rng(kSeed);
  std::bernoulli_distribution pick_two_qubit(0.5);
  std::uniform_int_distribution<std::size_t> one_qubit_gate(
      0, kOneQubit.size() - 1);
  std::uniform_int_distribution<std::size_t> two_qubit_gate(
      0, kTwoQubit.size() - 1);
  std::uniform_int_distribution<int> qubit(0, kQubits - 1);
  std::uniform_real_distribution<double> angle(0.0, 2 * std::numbers::pi);
  StateVector state(kQubits);
  for (int g = 0; g < kGates; ++g) {
    const bool two_qubit = pick_two_qubit(rng);
    const GateType gate = two_qubit ? kTwoQubit[two_qubit_gate(rng)]
                                    : kOneQubit[one_qubit_gate(rng)];
    const int j = qubit(rng);
    int k = j;  // second qubit; equal to j for a one-qubit gate
    if (two_qubit) {
      while (k == j) {
        k = qubit(rng);
      }
      state.apply_2q(matrix_2q(gate), j, k);
    } else {
      state.apply_1q(matrix_1q(gate, angle(rng)), j);
    }
    CAPTURE(g, static_cast<int>(gate), j, k);
    REQUIRE(state.probabilities().sum() == near(1.0));
  }
}

TEST_CASE("apply_2q rejects bad qubit pairs") {
  StateVector state(3);
  const Eigen::Matrix4cd cnot = matrix_2q(GateType::CNOT);
  REQUIRE_THROWS_AS(state.apply_2q(cnot, 1, 1), std::invalid_argument);
  REQUIRE_THROWS_AS(state.apply_2q(cnot, 0, 3), std::out_of_range);
  REQUIRE_THROWS_AS(state.apply_2q(cnot, -1, 0), std::out_of_range);
  REQUIRE_THROWS_AS(state.apply_2q(cnot, 3, 0), std::out_of_range);
}

TEST_CASE("apply_ry matches apply_1q with the RY matrix") {
  std::mt19937_64 rng(kSeed);
  for (const double theta : {0.0, 0.3, std::numbers::pi, 5.9}) {
    for (int n = 1; n <= 8; ++n) {
      for (int q = 0; q < n; ++q) {
        CAPTURE(theta, n, q);
        const Eigen::VectorXcd psi = random_state(n, rng).normalized();
        StateVector expected(psi);
        expected.apply_1q(matrix_1q(GateType::RY, theta), q);
        StateVector state(psi);
        state.apply_ry(theta, q);
        require_close(state.to_vector(), expected.to_vector());
      }
    }
  }
}

TEST_CASE("apply_cnot matches apply_2q with the CNOT matrix") {
  // Every ordered pair: control below and above the target, adjacent and
  // distant qubits.
  std::mt19937_64 rng(kSeed);
  for (int n = 2; n <= 8; ++n) {
    for (int c = 0; c < n; ++c) {
      for (int t = 0; t < n; ++t) {
        if (c == t) {
          continue;
        }
        CAPTURE(n, c, t);
        const Eigen::VectorXcd psi = random_state(n, rng).normalized();
        StateVector expected(psi);
        expected.apply_2q(matrix_2q(GateType::CNOT), c, t);
        StateVector state(psi);
        state.apply_cnot(c, t);
        require_close(state.to_vector(), expected.to_vector());
      }
    }
  }
}

TEST_CASE("apply_ry rejects a qubit outside the state") {
  const int q = GENERATE(-1, 3, 64);
  CAPTURE(q);
  StateVector state(3);
  REQUIRE_THROWS_AS(state.apply_ry(0.3, q), std::out_of_range);
}

TEST_CASE("apply_cnot rejects bad qubit pairs") {
  StateVector state(3);
  REQUIRE_THROWS_AS(state.apply_cnot(1, 1), std::invalid_argument);
  REQUIRE_THROWS_AS(state.apply_cnot(0, 3), std::out_of_range);
  REQUIRE_THROWS_AS(state.apply_cnot(-1, 0), std::out_of_range);
  REQUIRE_THROWS_AS(state.apply_cnot(3, 0), std::out_of_range);
  REQUIRE_THROWS_AS(state.apply_cnot(0, 64), std::out_of_range);
}

namespace {

// Real entries drawn from the standard normal distribution, as
// random_matrix draws its complex ones.
template <int Size>
Eigen::Matrix<double, Size, Size> random_real_matrix(std::mt19937_64& rng) {
  std::normal_distribution<double> normal;
  Eigen::Matrix<double, Size, Size> m;
  for (int r = 0; r < Size; ++r) {
    for (int c = 0; c < Size; ++c) {
      m(r, c) = normal(rng);
    }
  }
  return m;
}

}  // namespace

TEST_CASE("apply_real_1q matches apply_1q with the same matrix") {
  std::mt19937_64 rng(kSeed);
  for (int trial = 0; trial < 3; ++trial) {
    for (int n = 1; n <= 8; ++n) {
      for (int q = 0; q < n; ++q) {
        CAPTURE(trial, n, q);
        const Eigen::Matrix2d m = random_real_matrix<2>(rng);
        const Eigen::VectorXcd psi = random_state(n, rng);
        StateVector expected(psi);
        expected.apply_1q(m.cast<cd>(), q);
        StateVector state(psi);
        state.apply_real_1q(m, q);
        require_close(state.to_vector(), expected.to_vector());
      }
    }
  }
}

TEST_CASE("apply_real_2q matches apply_2q for every ordered qubit pair") {
  // first is the high bit of the matrix index in both kernels; pairs with
  // first below and above second cover both CNOT orientations.
  std::mt19937_64 rng(kSeed);
  for (int trial = 0; trial < 3; ++trial) {
    for (int n = 2; n <= 8; ++n) {
      for (int j = 0; j < n; ++j) {
        for (int k = 0; k < n; ++k) {
          if (j == k) {
            continue;
          }
          CAPTURE(trial, n, j, k);
          const Eigen::Matrix4d m = random_real_matrix<4>(rng);
          const Eigen::VectorXcd psi = random_state(n, rng);
          StateVector expected(psi);
          expected.apply_2q(m.cast<cd>(), j, k);
          StateVector state(psi);
          state.apply_real_2q(m, j, k);
          require_close(state.to_vector(), expected.to_vector());
        }
      }
    }
  }
}

TEST_CASE("apply_real_1q rejects a qubit outside the state") {
  const int k = GENERATE(-1, 3, 64);
  CAPTURE(k);
  StateVector state(3);
  const Eigen::Matrix2d x = matrix_1q(GateType::X).real();
  REQUIRE_THROWS_AS(state.apply_real_1q(x, k), std::out_of_range);
}

TEST_CASE("apply_real_2q rejects bad qubit pairs") {
  StateVector state(3);
  const Eigen::Matrix4d cnot = matrix_2q(GateType::CNOT).real();
  REQUIRE_THROWS_AS(state.apply_real_2q(cnot, 1, 1), std::invalid_argument);
  REQUIRE_THROWS_AS(state.apply_real_2q(cnot, 0, 3), std::out_of_range);
  REQUIRE_THROWS_AS(state.apply_real_2q(cnot, -1, 0), std::out_of_range);
  REQUIRE_THROWS_AS(state.apply_real_2q(cnot, 3, 0), std::out_of_range);
  REQUIRE_THROWS_AS(state.apply_real_2q(cnot, 0, 64), std::out_of_range);
}
