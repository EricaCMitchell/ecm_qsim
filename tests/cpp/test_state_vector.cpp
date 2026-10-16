#include <unistd.h>

#include <Eigen/Core>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <cmath>
#include <complex>
#include <cstdint>
#include <initializer_list>
#include <sstream>
#include <stdexcept>
#include <string>

#include "ecm_qsim/circuit.hpp"
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
  // 40 qubits need 2^40 * 24 bytes (24 TiB). For 63 qubits the byte count
  // does not fit in 64 bits.
  const int n = GENERATE(40, 63);
  CAPTURE(n);
  REQUIRE_THROWS_AS(StateVector(n), std::length_error);
  REQUIRE_THROWS_WITH(StateVector(n),
                      ContainsSubstring(std::to_string(n) + " qubits need 2^" +
                                        std::to_string(n) + " amplitudes"));
}

TEST_CASE("max_qubits fits 24 bytes per amplitude in half of memory") {
  // 2^limit amplitudes fit in the budget and 2^(limit + 1) do not. The
  // state at the limit is not built: it can take half of physical memory.
  const int limit = StateVector::max_qubits();
  CAPTURE(limit);
  REQUIRE(limit >= 20);  // 2^20 * 24 bytes is 24 MiB
  REQUIRE(limit < kMaxQubits);
  const std::uint64_t half_memory =
      static_cast<std::uint64_t>(sysconf(_SC_PHYS_PAGES)) *
      static_cast<std::uint64_t>(sysconf(_SC_PAGESIZE)) / 2;
  REQUIRE((std::uint64_t{24} << limit) <= half_memory);
  REQUIRE((std::uint64_t{24} << (limit + 1)) > half_memory);
  // The constructor refuses limit + 1 before it allocates.
  REQUIRE_THROWS_AS(StateVector(limit + 1), std::length_error);
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
