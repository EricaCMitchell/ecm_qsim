#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <map>
#include <numbers>
#include <stdexcept>
#include <string>
#include <vector>

#include "ecm_qsim/circuit.hpp"
#include "ecm_qsim/gates.hpp"
#include "ecm_qsim/measurement.hpp"
#include "ecm_qsim/noise.hpp"
#include "ecm_qsim/parser.hpp"
#include "ecm_qsim/simulator.hpp"
#include "ecm_qsim/state_vector.hpp"
#include "eigen_helpers.hpp"
#include "stat_helpers.hpp"

using namespace ecm_qsim;
using Catch::Matchers::ContainsSubstring;

namespace {

constexpr std::uint64_t kSeed = 42;

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

// The trajectory with one general-kernel call for each gate: apply_1q and
// apply_2q only, each step's errors after its gates.
StateVector general_kernels(const Circuit& circuit,
                            const std::vector<PauliError>& errors) {
  StateVector state(circuit.n_qubits);
  std::size_t next = 0;
  for (std::size_t t = 0; t < circuit.timesteps.size(); ++t) {
    for (const Op& op : circuit.timesteps[t].ops) {
      if (is_two_qubit(op.gate)) {
        state.apply_2q(matrix_2q(op.gate), op.qubits[0], op.qubits[1]);
      } else {
        state.apply_1q(matrix_1q(op.gate, op.theta), op.qubits[0]);
      }
    }
    for (; next < errors.size() && errors[next].timestep == t; ++next) {
      state.apply_1q(matrix_1q(errors[next].pauli), errors[next].qubit);
    }
  }
  return state;
}

// counts[bits], or 0 if bits never occurred.
int count_of(const Counts& counts, const std::string& bits) {
  const auto it = counts.find(bits);
  return it == counts.end() ? 0 : it->second;
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

TEST_CASE("simulate matches the general kernels on a noisy brickwork circuit") {
  // Three brickwork layers on 5 qubits (RY on every qubit, then CNOTs on
  // alternating neighboring pairs), plus a reversed CNOT and an H that take
  // the general paths, and X, Y and Z errors.
  constexpr int kQubits = 5;
  Circuit circuit{kQubits, {}};
  double theta = 0.1;
  for (int layer = 0; layer < 3; ++layer) {
    TimeStep rotations;
    for (int q = 0; q < kQubits; ++q) {
      rotations.ops.push_back(Op{GateType::RY, {q, 0}, theta});
      theta += 0.7;
    }
    TimeStep entanglers;
    for (int q = layer % 2; q + 1 < kQubits; q += 2) {
      entanglers.ops.push_back(Op{GateType::CNOT, {q, q + 1}, 0.0});
    }
    circuit.timesteps.push_back(rotations);
    circuit.timesteps.push_back(entanglers);
  }
  circuit.timesteps.push_back(TimeStep{
      {Op{GateType::CNOT, {4, 1}, 0.0}, Op{GateType::H, {2, 0}, 0.0}}});
  const std::vector<PauliError> errors{
      PauliError{0, 3, GateType::X}, PauliError{1, 0, GateType::Y},
      PauliError{1, 4, GateType::Z}, PauliError{4, 2, GateType::Y},
      PauliError{6, 1, GateType::X}};

  require_close(simulate(circuit, errors).to_vector(),
                general_kernels(circuit, errors).to_vector());
}

TEST_CASE("simulate matches the general kernels on depth-20 brickwork") {
  // Every RY and CNOT of the benchmark circuit, fused, against one
  // general-kernel call for each gate.
  constexpr int kQubits = 8;
  constexpr int kDepth = 20;
  Circuit circuit{kQubits, {}};
  double theta = 0.2;
  for (int layer = 1; layer <= kDepth; ++layer) {
    TimeStep rotations;
    for (int q = 0; q < kQubits; ++q) {
      rotations.ops.push_back(Op{GateType::RY, {q, 0}, theta});
      theta += 0.37;
    }
    TimeStep entanglers;
    for (int q = layer % 2 == 1 ? 0 : 1; q + 1 < kQubits; q += 2) {
      entanglers.ops.push_back(Op{GateType::CNOT, {q, q + 1}, 0.0});
    }
    circuit.timesteps.push_back(rotations);
    circuit.timesteps.push_back(entanglers);
  }
  const std::vector<std::vector<PauliError>> trajectories{
      {},
      {PauliError{0, 0, GateType::X}},
      {PauliError{3, 7, GateType::Y}, PauliError{3, 2, GateType::Z},
       PauliError{20, 4, GateType::X}},
      {PauliError{0, 1, GateType::Z}, PauliError{1, 1, GateType::Y},
       PauliError{17, 5, GateType::Y}, PauliError{39, 0, GateType::X},
       PauliError{39, 7, GateType::Z}},
  };
  for (std::size_t i = 0; i < trajectories.size(); ++i) {
    CAPTURE(i);
    require_close(simulate(circuit, trajectories[i]).to_vector(),
                  general_kernels(circuit, trajectories[i]).to_vector());
  }
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

TEST_CASE("simulate keeps its exception messages and their priority") {
  using Catch::Matchers::Equals;
  const auto one_op = [](int n, Op op) { return Circuit{n, {TimeStep{{op}}}}; };
  // Each gate reports a bad qubit with the name of its kernel.
  CHECK_THROWS_WITH(simulate(one_op(2, Op{GateType::RY, {2, 0}, 0.1}), {}),
                    Equals("apply_ry: qubit 2 out of range for 2 qubits"));
  CHECK_THROWS_WITH(simulate(one_op(2, Op{GateType::T, {5, 0}, 0.0}), {}),
                    Equals("apply_1q: qubit 5 out of range for 2 qubits"));
  CHECK_THROWS_WITH(simulate(one_op(2, Op{GateType::CNOT, {0, 2}, 0.0}), {}),
                    Equals("apply_cnot: qubit 2 out of range for 2 qubits"));
  CHECK_THROWS_WITH(simulate(one_op(2, Op{GateType::CNOT, {1, 1}, 0.0}), {}),
                    Equals("apply_cnot: both qubits are 1"));
  CHECK_THROWS_WITH(simulate(one_op(2, Op{GateType::CZ, {-1, 0}, 0.0}), {}),
                    Equals("apply_2q: qubit -1 out of range for 2 qubits"));
  CHECK_THROWS_WITH(simulate(one_op(2, Op{GateType::CY, {0, 0}, 0.0}), {}),
                    Equals("apply_2q: both qubits are 0"));
  CHECK_THROWS_WITH(
      simulate(Circuit{2, {TimeStep{}}}, {PauliError{0, 2, GateType::X}}),
      Equals("apply_1q: qubit 2 out of range for 2 qubits"));

  // Error metadata first, then the register size, then qubits in execution
  // order: a step's gates before its errors.
  const Circuit bad_gate = one_op(2, Op{GateType::RY, {2, 0}, 0.1});
  CHECK_THROWS_WITH(simulate(bad_gate, {PauliError{1, 0, GateType::X}}),
                    Equals("error timestep 1 out of range for 1 time steps"));
  CHECK_THROWS_AS(simulate(one_op(0, Op{GateType::RY, {2, 0}, 0.1}), {}),
                  std::invalid_argument);
  CHECK_THROWS_WITH(simulate(bad_gate, {PauliError{0, 7, GateType::X}}),
                    Equals("apply_ry: qubit 2 out of range for 2 qubits"));
  const Circuit late_gate{
      2, {TimeStep{}, TimeStep{{Op{GateType::RY, {2, 0}, 0.1}}}}};
  CHECK_THROWS_WITH(simulate(late_gate, {PauliError{0, 7, GateType::X}}),
                    Equals("apply_1q: qubit 7 out of range for 2 qubits"));
}

TEST_CASE("run rejects fewer than 1 shot") {
  const Circuit circuit{1, {}};
  const int bad = GENERATE(0, -1);
  CHECK_THROWS_AS(run(circuit, NoiseModel{}, bad, kSeed),
                  std::invalid_argument);
  CHECK_THROWS_WITH(run(circuit, NoiseModel{}, 0, kSeed),
                    "shots must be at least 1, got 0");
}

TEST_CASE("run rejects a bad noise model") {
  CHECK_THROWS_WITH(run(Circuit{1, {}}, NoiseModel{1.5, 0, 0}, 1, kSeed),
                    ContainsSubstring("p_gate must be in [0, 1]"));
}

TEST_CASE("run reports qubit 0 as the rightmost character") {
  const Circuit circuit{3, {TimeStep{{Op{GateType::X, {0, 0}, 0.0}}}}};
  CHECK(run(circuit, NoiseModel{}, 100, kSeed) == Counts{{"001", 100}});
}

TEST_CASE("run without noise gives the ideal GHZ distribution") {
  const Circuit circuit = parse_circuit_file(ECM_QSIM_EXAMPLES_DIR "/ghz3.txt");
  const Counts counts = run(circuit, NoiseModel{}, 10000, 1);
  REQUIRE(counts.size() == 2);
  REQUIRE(count_of(counts, "000") + count_of(counts, "111") == 10000);
  require_binomial(count_of(counts, "000"), 10000, 0.5);
}

TEST_CASE("run with readout error only returns 0 after X at rate p_readout") {
  const Circuit circuit{1, {TimeStep{{Op{GateType::X, {0, 0}, 0.0}}}}};
  const Counts counts = run(circuit, NoiseModel{0, 0, 0.1}, 100000, kSeed);
  REQUIRE(count_of(counts, "0") + count_of(counts, "1") == 100000);
  require_binomial(count_of(counts, "0"), 100000, 0.1);
}

TEST_CASE(
    "run gives gate noise to the gate's qubit and idle noise to the others") {
  // X 0 on 2 qubits. An X or Y error flips the measured bit (rate 2/3); a Z
  // error does not. A state kept from the previous shot would give about 1/2.
  constexpr int kShots = 30000;
  const Circuit circuit{2, {TimeStep{{Op{GateType::X, {0, 0}, 0.0}}}}};
  SECTION("p_gate = 1 acts on qubit 0 only") {
    const Counts counts = run(circuit, NoiseModel{1, 0, 0}, kShots, kSeed);
    for (const auto& [bits, n] : counts) {
      CAPTURE(bits);
      CHECK((bits == "00" || bits == "01"));
    }
    require_binomial(count_of(counts, "00"), kShots, 2.0 / 3);
  }
  SECTION("p_idle = 1 acts on qubit 1 only") {
    const Counts counts = run(circuit, NoiseModel{0, 1, 0}, kShots, kSeed);
    for (const auto& [bits, n] : counts) {
      CAPTURE(bits);
      CHECK((bits == "01" || bits == "11"));
    }
    require_binomial(count_of(counts, "11"), kShots, 2.0 / 3);
  }
}

// This pins the Part 1 reference baseline only. Optimised engines may reorder
// draws and will be checked for statistical equivalence instead.
TEST_CASE("run reference engine follows the documented draw order") {
  const Circuit circuit = parse_circuit_file(ECM_QSIM_EXAMPLES_DIR "/ghz3.txt");
  const NoiseModel noise{0.05, 0.05, 0.05};
  constexpr int kShots = 500;
  Rng rng(7);
  Counts expected;
  for (int shot = 0; shot < kShots; ++shot) {
    const auto errors = sample_pauli_errors(circuit, noise, rng);
    const std::uint64_t index = sample_index(
        simulate(circuit, errors).probabilities(), uniform_from_bits(rng()));
    ++expected[to_bitstring(
        flip_readout_bits(index, circuit.n_qubits, noise.p_readout, rng),
        circuit.n_qubits)];
  }
  REQUIRE(expected.size() > 2);  // the noise did fire
  CHECK(run(circuit, noise, kShots, 7) == expected);
}

TEST_CASE("run is reproducible from the seed") {
  const Circuit circuit = parse_circuit_file(ECM_QSIM_EXAMPLES_DIR "/ghz3.txt");
  const NoiseModel noise{0.05, 0.05, 0.05};
  const Counts first = run(circuit, noise, 1000, 1);
  CHECK(run(circuit, noise, 1000, 1) == first);
  CHECK(run(circuit, noise, 1000, 2) != first);
}

TEST_CASE("run refuses a register that does not fit in memory") {
  const int n = StateVector::max_qubits() + 1;
  if (n > kMaxQubits) {
    SKIP("this machine fits every register the IR can represent");
  }
  CHECK_THROWS_AS(run(Circuit{n, {}}, NoiseModel{}, 1, kSeed),
                  std::length_error);
}
