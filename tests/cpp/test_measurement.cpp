#include <Eigen/Core>
#include <array>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <cmath>
#include <complex>
#include <cstdint>
#include <initializer_list>
#include <limits>
#include <map>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include "ecm_qsim/circuit.hpp"
#include "ecm_qsim/gates.hpp"
#include "ecm_qsim/measurement.hpp"
#include "ecm_qsim/state_vector.hpp"

using namespace ecm_qsim;

namespace {

Eigen::VectorXd probs(std::initializer_list<double> entries) {
  Eigen::VectorXd p(static_cast<Eigen::Index>(entries.size()));
  Eigen::Index i = 0;
  for (const double x : entries) {
    p(i++) = x;
  }
  return p;
}

// A normalized state on n_qubits qubits with deterministic random amplitudes.
StateVector random_state(int n_qubits, std::uint64_t seed) {
  std::mt19937_64 rng(seed);
  std::normal_distribution<double> normal;
  Eigen::VectorXcd amplitudes(Eigen::Index{1} << n_qubits);
  for (Eigen::Index i = 0; i < amplitudes.size(); ++i) {
    amplitudes(i) = {normal(rng), normal(rng)};
  }
  amplitudes.normalize();
  return StateVector(amplitudes);
}

// The what() of the std::invalid_argument that call throws.
template <typename Call>
std::string invalid_argument_message(Call call) {
  try {
    call();
  } catch (const std::invalid_argument& e) {
    return e.what();
  }
  FAIL("no std::invalid_argument was thrown");
  return {};
}

}  // namespace

TEST_CASE("sample_index maps u onto the cumulative distribution") {
  const Eigen::VectorXd p = probs({0, 0.25, 0, 0.75});
  REQUIRE(sample_index(p, 0.0) == 1);
  REQUIRE(sample_index(p, 0.2) == 1);
  REQUIRE(sample_index(p, 0.3) == 3);
  REQUIRE(sample_index(p, std::nextafter(1.0, 0.0)) == 3);
}

TEST_CASE("sample_index never returns a zero-probability index") {
  // Zeros at both ends: u = 0 must skip index 0, and u just below 1 must
  // not run past index 2.
  const Eigen::VectorXd p = probs({0, 0.36, 0.64, 0});
  REQUIRE(sample_index(p, 0.0) == 1);
  REQUIRE(sample_index(p, std::nextafter(1.0, 0.0)) == 2);
}

TEST_CASE("sample_index puts a boundary value in the next interval") {
  // Index i owns [c(i - 1), c(i)), where c is the running sum. All values
  // here are exact in binary.
  const Eigen::VectorXd p = probs({0.25, 0.25, 0.5});
  REQUIRE(sample_index(p, 0.25) == 1);
  REQUIRE(sample_index(p, 0.5) == 2);
}

TEST_CASE("sample_index falls back to the last non-zero index") {
  // With a subnormal total, u * total rounds up to total for u just below
  // 1, so no running sum is greater. The result must still have p > 0.
  const double tiny = std::numeric_limits<double>::denorm_min();
  const Eigen::VectorXd p = probs({0, tiny, 0});
  REQUIRE(sample_index(p, std::nextafter(1.0, 0.0)) == 1);
}

TEST_CASE("sample_index scales u by the total probability") {
  const Eigen::VectorXd p = probs({1, 1});  // not normalised: total is 2
  REQUIRE(sample_index(p, 0.49) == 0);
  REQUIRE(sample_index(p, 0.51) == 1);
}

TEST_CASE("sample_index reproduces the distribution over an even grid of u") {
  // Grid points sit 0.0005 away from each boundary.
  const Eigen::VectorXd p = probs({0.2, 0.3, 0.5, 0});
  constexpr int kDraws = 1000;
  std::array<int, 4> counts{};
  for (int t = 0; t < kDraws; ++t) {
    ++counts[sample_index(p, (t + 0.5) / kDraws)];
  }
  REQUIRE(counts == std::array<int, 4>{200, 300, 500, 0});
}

TEST_CASE("GHZ samples give only 000 and 111, half each") {
  StateVector state(3);
  state.apply_1q(matrix_1q(GateType::H), 0);
  state.apply_2q(matrix_2q(GateType::CNOT), 0, 1);
  state.apply_2q(matrix_2q(GateType::CNOT), 1, 2);
  const Eigen::VectorXd p = state.probabilities();
  constexpr int kDraws = 1000;
  std::map<std::string, int> counts;
  for (int t = 0; t < kDraws; ++t) {
    ++counts[to_bitstring(sample_index(p, (t + 0.5) / kDraws), 3)];
  }
  REQUIRE(counts == std::map<std::string, int>{{"000", 500}, {"111", 500}});
}

TEST_CASE("sample_index rejects u outside 0 to 1, and NaN") {
  const double u =
      GENERATE(-0.1, 1.0, 1.5, std::numeric_limits<double>::quiet_NaN());
  CAPTURE(u);
  REQUIRE_THROWS_AS(sample_index(probs({1, 0}), u), std::invalid_argument);
}

TEST_CASE("sample_index rejects negative, NaN and infinite entries") {
  const double bad = GENERATE(-0.1, std::numeric_limits<double>::quiet_NaN(),
                              std::numeric_limits<double>::infinity());
  CAPTURE(bad);
  REQUIRE_THROWS_AS(sample_index(probs({0.5, bad}), 0.5),
                    std::invalid_argument);
}

TEST_CASE("sample_index rejects a distribution with no weight") {
  REQUIRE_THROWS_AS(sample_index(probs({0, 0}), 0.5), std::domain_error);
  REQUIRE_THROWS_AS(sample_index(Eigen::VectorXd(), 0.5), std::domain_error);
}

TEST_CASE("sample_index on a state agrees with its probabilities") {
  const int n_qubits = GENERATE(1, 6);
  const std::uint64_t seed = GENERATE(1, 2, 3, 4);
  CAPTURE(n_qubits, seed);
  const StateVector state = random_state(n_qubits, seed);
  const Eigen::VectorXd p = state.probabilities();

  // Values of u just below, at and just above each cumulative boundary,
  // summed in the order sample_index adds them, plus an even grid.
  const double total = p.sum();
  std::vector<double> us;
  double cumulative = 0.0;
  for (Eigen::Index i = 0; i < p.size(); ++i) {
    cumulative += p(i);
    const double b = cumulative / total;
    for (const double u : {std::nextafter(b, 0.0), b, std::nextafter(b, 1.0)}) {
      if (u >= 0.0 && u < 1.0) {
        us.push_back(u);
      }
    }
  }
  constexpr int kGrid = 1024;
  for (int t = 0; t < kGrid; ++t) {
    us.push_back(static_cast<double>(t) / kGrid);
  }

  for (const double u : us) {
    CAPTURE(u);
    REQUIRE(sample_index(state, u) == sample_index(p, u));
  }
}

TEST_CASE("sample_index on a state handles exceptional amplitudes") {
  SECTION("a subnormal probability reaches the last non-zero fallback") {
    Eigen::VectorXcd amplitudes = Eigen::VectorXcd::Zero(4);
    amplitudes(1) = 1e-160;  // |amplitude|^2 = 1e-320 is subnormal
    REQUIRE(sample_index(StateVector(amplitudes), std::nextafter(1.0, 0.0)) ==
            1);
  }
  SECTION("NaN and infinite amplitudes throw the vector overload's error") {
    const double bad = GENERATE(std::numeric_limits<double>::quiet_NaN(),
                                std::numeric_limits<double>::infinity());
    CAPTURE(bad);
    Eigen::VectorXcd amplitudes(2);
    amplitudes << 0.5, bad;
    const StateVector state(amplitudes);
    const std::string expected = invalid_argument_message(
        [&] { sample_index(state.probabilities(), 0.5); });
    REQUIRE(invalid_argument_message([&] { sample_index(state, 0.5); }) ==
            expected);
    REQUIRE(expected ==
            (std::isnan(bad)
                 ? "sample_index: probability 1 is nan"
                 : "sample_index: total probability is not finite"));
  }
}

TEST_CASE("uniform_from_bits maps 64 bits into [0, 1)") {
  REQUIRE(uniform_from_bits(0) == 0.0);
  REQUIRE(uniform_from_bits(std::uint64_t{1} << 63) == 0.5);
  // The largest input gives the largest double below 1, never 1.0.
  REQUIRE(uniform_from_bits(std::numeric_limits<std::uint64_t>::max()) ==
          std::nextafter(1.0, 0.0));
}
