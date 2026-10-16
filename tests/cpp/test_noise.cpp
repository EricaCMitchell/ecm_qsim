#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <cmath>
#include <cstdint>
#include <initializer_list>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include "ecm_qsim/circuit.hpp"
#include "ecm_qsim/noise.hpp"

using namespace ecm_qsim;
using Catch::Matchers::ContainsSubstring;

namespace {

constexpr std::uint64_t kSeed = 42;
constexpr int kTrials = 100000;

// |count - trials * p| within 5 binomial standard deviations.
void require_binomial(int count, int trials, double p) {
  const double mean = trials * p;
  const double sigma = std::sqrt(trials * p * (1 - p));
  CAPTURE(count, trials, p, mean, sigma);
  REQUIRE(std::abs(count - mean) <= 5 * sigma);
}

// In both circuits the gate's qubits are busy and the other qubits are idle.
// kH2Circuit leaves qubit 0 idle next to a one-qubit gate, which catches a
// sampler that marks the unused qubits[1] = 0.
const Circuit kCnotCircuit{3, {TimeStep{{Op{GateType::CNOT, {0, 1}, 0.0}}}}};
const Circuit kH2Circuit{3, {TimeStep{{Op{GateType::H, {2, 0}, 0.0}}}}};

}  // namespace

namespace Catch {
template <>
struct StringMaker<ecm_qsim::PauliError> {
  static std::string convert(const ecm_qsim::PauliError& error) {
    return "PauliError{timestep=" + std::to_string(error.timestep) +
           ", qubit=" + std::to_string(error.qubit) +
           ", pauli=" + std::to_string(static_cast<int>(error.pauli)) + "}";
  }
};
}  // namespace Catch

TEST_CASE("NoiseModel defaults to noiseless") {
  NoiseModel noise;
  REQUIRE(noise.p_gate == 0.0);
  REQUIRE(noise.p_idle == 0.0);
  REQUIRE(noise.p_readout == 0.0);
}

TEST_CASE("validate_noise_model accepts probabilities in [0, 1]") {
  CHECK_NOTHROW(validate_noise_model(NoiseModel{}));
  CHECK_NOTHROW(validate_noise_model(NoiseModel{1.0, 1.0, 1.0}));
  CHECK_NOTHROW(validate_noise_model(NoiseModel{0.001, 0.5, 0.01}));
}

TEST_CASE("validate_noise_model rejects values outside [0, 1] and NaN") {
  const double nan = std::numeric_limits<double>::quiet_NaN();
  const double inf = std::numeric_limits<double>::infinity();
  const double bad = GENERATE_COPY(-0.1, 1.5, nan, inf, -inf);
  CAPTURE(bad);
  CHECK_THROWS_WITH(validate_noise_model(NoiseModel{bad, 0, 0}),
                    ContainsSubstring("p_gate must be in [0, 1], got "));
  CHECK_THROWS_WITH(validate_noise_model(NoiseModel{0, bad, 0}),
                    ContainsSubstring("p_idle must be in [0, 1], got "));
  CHECK_THROWS_WITH(validate_noise_model(NoiseModel{0, 0, bad}),
                    ContainsSubstring("p_readout must be in [0, 1], got "));
  CHECK_THROWS_AS(validate_noise_model(NoiseModel{bad, 0, 0}),
                  std::invalid_argument);
}

TEST_CASE("validate_noise_model reports the value as given") {
  CHECK_THROWS_WITH(validate_noise_model(NoiseModel{0, 1.5, 0}),
                    "p_idle must be in [0, 1], got 1.5");
  // Default ostream formatting would print "got 1" here.
  CHECK_THROWS_WITH(validate_noise_model(NoiseModel{1.0000001, 0, 0}),
                    "p_gate must be in [0, 1], got 1.0000001");
  CHECK_THROWS_WITH(validate_noise_model(NoiseModel{
                        0, 0, -std::numeric_limits<double>::infinity()}),
                    "p_readout must be in [0, 1], got -inf");
}

TEST_CASE("validate_noise_model reports the first bad field") {
  const double nan = std::numeric_limits<double>::quiet_NaN();
  CHECK_THROWS_WITH(validate_noise_model(NoiseModel{-0.1, 1.5, nan}),
                    "p_gate must be in [0, 1], got -0.1");
  CHECK_THROWS_WITH(validate_noise_model(NoiseModel{0, 1.5, nan}),
                    "p_idle must be in [0, 1], got 1.5");
}

TEST_CASE("sample_pauli_errors gives no errors and no draws when noiseless") {
  Rng rng(kSeed);
  const Rng before = rng;
  CHECK(sample_pauli_errors(kCnotCircuit, NoiseModel{0, 0, 0.5}, rng).empty());
  CHECK(rng == before);
}

TEST_CASE("Gate noise hits both qubits of a two-qubit gate and no idle qubit") {
  Rng rng(kSeed);
  const auto errors =
      sample_pauli_errors(kCnotCircuit, NoiseModel{1, 0, 0}, rng);
  REQUIRE(errors.size() == 2);
  CHECK(errors[0].timestep == 0);
  CHECK(errors[0].qubit == 0);
  CHECK(errors[1].timestep == 0);
  CHECK(errors[1].qubit == 1);
}

TEST_CASE("Idle noise hits only qubits with no gate in the step") {
  Rng rng(kSeed);
  const auto errors =
      sample_pauli_errors(kCnotCircuit, NoiseModel{0, 1, 0}, rng);
  REQUIRE(errors.size() == 1);
  CHECK(errors[0].qubit == 2);
}

TEST_CASE("Gate noise on a one-qubit gate ignores the unused second qubit") {
  Rng rng(kSeed);
  const auto errors = sample_pauli_errors(kH2Circuit, NoiseModel{1, 0, 0}, rng);
  REQUIRE(errors.size() == 1);
  CHECK(errors[0].qubit == 2);
}

TEST_CASE("Idle noise hits qubit 0 when a one-qubit gate acts elsewhere") {
  Rng rng(kSeed);
  const auto errors = sample_pauli_errors(kH2Circuit, NoiseModel{0, 1, 0}, rng);
  REQUIRE(errors.size() == 2);
  CHECK(errors[0].qubit == 0);
  CHECK(errors[1].qubit == 1);
}

TEST_CASE("Every qubit is idle in an empty step; errors are sorted") {
  const Circuit circuit{2, {TimeStep{}, TimeStep{}, TimeStep{}}};
  Rng rng(kSeed);
  const auto errors = sample_pauli_errors(circuit, NoiseModel{0, 1, 0}, rng);
  REQUIRE(errors.size() == 6);
  for (std::size_t i = 0; i < errors.size(); ++i) {
    CAPTURE(i);
    CHECK(errors[i].timestep == i / 2);
    CHECK(errors[i].qubit == static_cast<int>(i % 2));
    CHECK((errors[i].pauli == GateType::X || errors[i].pauli == GateType::Y ||
           errors[i].pauli == GateType::Z));
  }
}

TEST_CASE("A circuit with no time steps has no errors") {
  Rng rng(kSeed);
  CHECK(sample_pauli_errors(Circuit{4, {}}, NoiseModel{1, 1, 1}, rng).empty());
}

TEST_CASE("sample_pauli_errors rejects invalid circuit qubit counts") {
  Rng rng(kSeed);
  const int bad = GENERATE(-1, 0, kMaxQubits + 1);
  CAPTURE(bad);
  CHECK_THROWS_AS(sample_pauli_errors(Circuit{bad, {}}, NoiseModel{}, rng),
                  std::invalid_argument);
}

TEST_CASE("sample_pauli_errors rejects a bad qubit before any draw") {
  const int bad = GENERATE(-1, 2);
  CAPTURE(bad);
  // The first step is valid and noisy, so a check made during sampling
  // would draw before it throws.
  const Circuit bad_1q{
      2, {TimeStep{}, TimeStep{{Op{GateType::H, {bad, 0}, 0.0}}}}};
  const Circuit bad_2q{
      2, {TimeStep{}, TimeStep{{Op{GateType::CNOT, {0, bad}, 0.0}}}}};
  for (const Circuit& circuit : {bad_1q, bad_2q}) {
    Rng rng(kSeed);
    const Rng before = rng;
    CHECK_THROWS_AS(sample_pauli_errors(circuit, NoiseModel{1, 1, 0}, rng),
                    std::out_of_range);
    CHECK(rng == before);
  }
}

TEST_CASE("sample_pauli_errors does not check qubits[1] of a one-qubit gate") {
  Rng rng(kSeed);
  const Circuit circuit{2, {TimeStep{{Op{GateType::H, {0, 7}, 0.0}}}}};
  CHECK_NOTHROW(sample_pauli_errors(circuit, NoiseModel{1, 1, 0}, rng));
}

TEST_CASE("Gate and idle errors occur at their own rates") {
  const Circuit circuit{2, {TimeStep{{Op{GateType::H, {0, 0}, 0.0}}}}};
  Rng rng(kSeed);
  int gate_count = 0;
  int idle_count = 0;
  for (int trial = 0; trial < kTrials; ++trial) {
    for (const PauliError& error :
         sample_pauli_errors(circuit, NoiseModel{0.1, 0.4, 0}, rng)) {
      if (error.qubit == 0) {
        ++gate_count;
      } else {
        ++idle_count;
      }
    }
  }
  require_binomial(gate_count, kTrials, 0.1);
  require_binomial(idle_count, kTrials, 0.4);
}

TEST_CASE("Pauli errors are X, Y and Z with equal probability") {
  const Circuit circuit{1, {TimeStep{}}};
  Rng rng(kSeed);
  int x_count = 0;
  int y_count = 0;
  int z_count = 0;
  for (int trial = 0; trial < kTrials; ++trial) {
    const auto errors = sample_pauli_errors(circuit, NoiseModel{0, 1, 0}, rng);
    REQUIRE(errors.size() == 1);
    switch (errors[0].pauli) {
      case GateType::X:
        ++x_count;
        break;
      case GateType::Y:
        ++y_count;
        break;
      case GateType::Z:
        ++z_count;
        break;
      default:
        FAIL("pauli is not X, Y or Z");
    }
  }
  require_binomial(x_count, kTrials, 1.0 / 3);
  require_binomial(y_count, kTrials, 1.0 / 3);
  require_binomial(z_count, kTrials, 1.0 / 3);
}

TEST_CASE("The two gate-error draws of a two-qubit gate are independent") {
  // p_gate = 0.5: both qubits hit at rate 1/4, and both hit with the same
  // Pauli at rate 1/12. One shared occurrence draw gives 1/2; one shared
  // Pauli gives 1/4.
  const Circuit circuit{2, {TimeStep{{Op{GateType::CNOT, {0, 1}, 0.0}}}}};
  Rng rng(kSeed);
  int both = 0;
  int same_pauli = 0;
  for (int trial = 0; trial < kTrials; ++trial) {
    const auto errors =
        sample_pauli_errors(circuit, NoiseModel{0.5, 0, 0}, rng);
    if (errors.size() == 2) {
      ++both;
      if (errors[0].pauli == errors[1].pauli) {
        ++same_pauli;
      }
    }
  }
  require_binomial(both, kTrials, 0.25);
  require_binomial(same_pauli, kTrials, 0.25 / 3);
}

TEST_CASE("sample_pauli_errors uses the fixed draw order") {
  // Step 0: CNOT 0 1, qubit 2 idle. Step 1: H 2, qubits 0 and 1 idle.
  // Step 2: empty. 9 occurrence draws + 5 Pauli draws = 14. The CNOT's
  // control and target get different Paulis (Z, X).
  const Circuit circuit{3,
                        {TimeStep{{Op{GateType::CNOT, {0, 1}, 0.0}}},
                         TimeStep{{Op{GateType::H, {2, 0}, 0.0}}}, TimeStep{}}};
  Rng rng(kSeed);
  const auto errors =
      sample_pauli_errors(circuit, NoiseModel{0.8, 0.4, 0}, rng);
  const std::vector<PauliError> expected_errors{{0, 0, GateType::Z},
                                                {0, 1, GateType::X},
                                                {1, 0, GateType::Y},
                                                {1, 1, GateType::Y},
                                                {1, 2, GateType::Y}};
  CHECK(errors == expected_errors);

  Rng expected_rng(kSeed);
  for (int draw = 0; draw < 14; ++draw) {
    static_cast<void>(expected_rng());
  }
  CHECK(rng == expected_rng);
}

TEST_CASE("sample_pauli_errors rejects an invalid noise model") {
  Rng rng(kSeed);
  CHECK_THROWS_WITH(
      sample_pauli_errors(kCnotCircuit, NoiseModel{0, std::nan(""), 0}, rng),
      ContainsSubstring("p_idle must be in [0, 1]"));
}

TEST_CASE("flip_readout_bits with p = 0 returns the index and draws nothing") {
  Rng rng(kSeed);
  const Rng before = rng;
  CHECK(flip_readout_bits(0b101, 3, 0.0, rng) == 0b101);
  CHECK(rng == before);
}

TEST_CASE("flip_readout_bits with p = 1 flips every bit") {
  Rng rng(kSeed);
  CHECK(flip_readout_bits(0b101, 3, 1.0, rng) == 0b010);
  CHECK(flip_readout_bits(0, 1, 1.0, rng) == 1);
  CHECK(flip_readout_bits(0, kMaxQubits, 1.0, rng) ==
        (std::uint64_t{1} << 63) - 1);
}

TEST_CASE("flip_readout_bits flips each bit at rate p") {
  constexpr std::uint64_t kIndex = 0b101;
  Rng rng(kSeed);
  int counts[3] = {0, 0, 0};
  for (int trial = 0; trial < kTrials; ++trial) {
    const std::uint64_t flipped =
        flip_readout_bits(kIndex, 3, 0.25, rng) ^ kIndex;
    for (int k = 0; k < 3; ++k) {
      if ((flipped >> k) & 1) {
        ++counts[k];
      }
    }
  }
  for (int k = 0; k < 3; ++k) {
    CAPTURE(k);
    require_binomial(counts[k], kTrials, 0.25);
  }
}

TEST_CASE("flip_readout_bits flips bits independently") {
  // One draw that flips all bits together gives 0.25.
  constexpr std::uint64_t kIndex = 0b101;
  Rng rng(kSeed);
  int both = 0;
  for (int trial = 0; trial < kTrials; ++trial) {
    const std::uint64_t flipped =
        flip_readout_bits(kIndex, 3, 0.25, rng) ^ kIndex;
    if ((flipped & 0b101) == 0b101) {
      ++both;
    }
  }
  require_binomial(both, kTrials, 0.25 * 0.25);
}

TEST_CASE("flip_readout_bits uses the fixed draw order") {
  // 5 draws; only the 4th (bit 3) fires. A reversed draw order flips bit 1.
  Rng rng(kSeed);
  CHECK(flip_readout_bits(0b10110, 5, 0.5, rng) == 0b11110);
  Rng expected_rng(kSeed);
  for (int draw = 0; draw < 5; ++draw) {
    static_cast<void>(expected_rng());
  }
  CHECK(rng == expected_rng);
}

TEST_CASE("flip_readout_bits rejects bad arguments") {
  Rng rng(kSeed);
  CHECK_THROWS_WITH(flip_readout_bits(0, 3, 1.5, rng),
                    "p_readout must be in [0, 1], got 1.5");
  CHECK_THROWS_AS(flip_readout_bits(0, 3, std::nan(""), rng),
                  std::invalid_argument);
  CHECK_THROWS_AS(flip_readout_bits(0, 0, 0.1, rng), std::invalid_argument);
  CHECK_THROWS_AS(flip_readout_bits(0, kMaxQubits + 1, 0.1, rng),
                  std::invalid_argument);
  CHECK_THROWS_AS(flip_readout_bits(8, 3, 0.1, rng), std::out_of_range);
}
