#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <cmath>
#include <complex>
#include <cstdint>
#include <numbers>
#include <stdexcept>
#include <vector>

#include "ecm_qsim/circuit.hpp"
#include "ecm_qsim/noise.hpp"
#include "ecm_qsim/parser.hpp"
#include "ecm_qsim/simulator.hpp"
#include "ecm_qsim/state_vector.hpp"
#include "eigen_helpers.hpp"

using namespace ecm_qsim;
using Catch::Matchers::ContainsSubstring;

namespace {

// |amplitude(index)| is 1 and every other amplitude is 0 (within 1e-12),
// and amplitude(index) equals expected.
void require_basis_state(const StateVector& state, std::uint64_t index,
                         std::complex<double> expected = 1.0) {
  for (std::uint64_t i = 0; i < state.size(); ++i) {
    CAPTURE(i);
    const std::complex<double> want = i == index ? expected : 0.0;
    REQUIRE(state.amplitude(i).real() == near(want.real()));
    REQUIRE(state.amplitude(i).imag() == near(want.imag()));
  }
}

}  // namespace

TEST_CASE("simulate of a circuit with no time steps gives |0...0>") {
  const StateVector state = simulate(Circuit{3, {}}, {});
  REQUIRE(state.n_qubits() == 3);
  require_basis_state(state, 0);
}

TEST_CASE("simulate applies a one-qubit gate to qubits[0] only") {
  // qubits[1] is 0 and unused; touching it would also act on qubit 0.
  const Circuit circuit{3, {TimeStep{{Op{GateType::X, {2, 0}, 0.0}}}}};
  require_basis_state(simulate(circuit, {}), 0b100);
}

TEST_CASE("simulate passes the angle of rotation gates") {
  const Circuit circuit{
      1, {TimeStep{{Op{GateType::RY, {0, 0}, std::numbers::pi}}}}};
  require_basis_state(simulate(circuit, {}), 1);  // RY(pi)|0> = |1>
}

TEST_CASE("simulate passes the control as the first qubit") {
  const TimeStep x0{{Op{GateType::X, {0, 0}, 0.0}}};
  SECTION("CNOT 0 1 flips qubit 1") {
    const Circuit c{2, {x0, TimeStep{{Op{GateType::CNOT, {0, 1}, 0.0}}}}};
    require_basis_state(simulate(c, {}), 0b11);
  }
  SECTION("CNOT 1 0 does nothing when qubit 1 is 0") {
    const Circuit c{2, {x0, TimeStep{{Op{GateType::CNOT, {1, 0}, 0.0}}}}};
    require_basis_state(simulate(c, {}), 0b01);
  }
}

TEST_CASE("simulate runs the example GHZ circuit, including the last step") {
  const Circuit circuit = parse_circuit_file(ECM_QSIM_EXAMPLES_DIR "/ghz3.txt");
  const StateVector state = simulate(circuit, {});
  // RZ(0.785) on qubit 0 in the final step, which has no trailing TICK.
  const double half = 0.785 / 2;
  const std::complex<double> a0 = std::polar(1 / std::sqrt(2.0), -half);
  const std::complex<double> a7 = std::polar(1 / std::sqrt(2.0), half);
  for (std::uint64_t i = 0; i < 8; ++i) {
    CAPTURE(i);
    const std::complex<double> expected = i == 0   ? a0
                                          : i == 7 ? a7
                                                   : std::complex<double>{};
    CHECK(state.amplitude(i).real() == near(expected.real()));
    CHECK(state.amplitude(i).imag() == near(expected.imag()));
  }
}

TEST_CASE("simulate applies an error after the gates of its step") {
  // H Z H |0> = |1>. Z before the first H, or after the second, gives |0>.
  const TimeStep h0{{Op{GateType::H, {0, 0}, 0.0}}};
  const Circuit circuit{1, {h0, h0}};
  require_basis_state(simulate(circuit, {}), 0);
  require_basis_state(simulate(circuit, {PauliError{0, 0, GateType::Z}}), 1);
}

TEST_CASE("simulate applies errors in an empty time step") {
  const Circuit circuit{2, {TimeStep{}}};
  require_basis_state(simulate(circuit, {PauliError{0, 1, GateType::X}}), 0b10);
}

TEST_CASE("simulate applies the Pauli matrix of each error, in step order") {
  const Circuit circuit{1, {TimeStep{}, TimeStep{}}};
  // Y|0> = i|1>.
  require_basis_state(simulate(circuit, {PauliError{0, 0, GateType::Y}}), 1,
                      {0.0, 1.0});
  // Z X |0> = -|1>; the reverse order X Z |0> gives +|1>.
  require_basis_state(simulate(circuit, {PauliError{0, 0, GateType::X},
                                         PauliError{1, 0, GateType::Z}}),
                      1, -1.0);
}

TEST_CASE("simulate rejects bad errors") {
  const Circuit circuit{3, {TimeStep{}, TimeStep{}}};
  const std::vector<PauliError> late{PauliError{2, 0, GateType::X}};
  CHECK_THROWS_AS(simulate(circuit, late), std::out_of_range);
  CHECK_THROWS_WITH(
      simulate(circuit, late),
      ContainsSubstring("error timestep 2 out of range for 2 time steps"));

  const std::vector<PauliError> unsorted{PauliError{1, 0, GateType::X},
                                         PauliError{0, 0, GateType::X}};
  CHECK_THROWS_AS(simulate(circuit, unsorted), std::invalid_argument);
  CHECK_THROWS_WITH(simulate(circuit, unsorted),
                    ContainsSubstring("errors must be sorted by timestep"));

  const std::vector<PauliError> not_pauli{PauliError{0, 0, GateType::H}};
  CHECK_THROWS_AS(simulate(circuit, not_pauli), std::invalid_argument);
  CHECK_THROWS_WITH(simulate(circuit, not_pauli),
                    ContainsSubstring("error pauli must be X, Y or Z"));

  CHECK_THROWS_AS(simulate(circuit, {PauliError{0, 3, GateType::X}}),
                  std::out_of_range);
}

TEST_CASE("simulate rejects bad circuits") {
  CHECK_THROWS_AS(simulate(Circuit{0, {}}, {}), std::invalid_argument);
  const Circuit bad_qubit{2, {TimeStep{{Op{GateType::X, {2, 0}, 0.0}}}}};
  CHECK_THROWS_AS(simulate(bad_qubit, {}), std::out_of_range);
}
