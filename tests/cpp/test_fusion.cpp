#include <Eigen/Core>
#include <array>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <numbers>
#include <random>
#include <stdexcept>
#include <variant>
#include <vector>

#include "ecm_qsim/circuit.hpp"
#include "ecm_qsim/fusion.hpp"
#include "ecm_qsim/gates.hpp"
#include "ecm_qsim/noise.hpp"
#include "ecm_qsim/state_vector.hpp"
#include "eigen_helpers.hpp"

using namespace ecm_qsim;
using Catch::Matchers::Equals;

namespace {

using cd = std::complex<double>;

constexpr std::uint64_t kSeed = 42;

constexpr std::array<GateType, 12> kAllGates = {
    GateType::X,  GateType::Y,    GateType::Z,  GateType::H,
    GateType::S,  GateType::T,    GateType::RX, GateType::RY,
    GateType::RZ, GateType::CNOT, GateType::CZ, GateType::CY};

// The trajectory with only the general kernels: for each step, its gates in
// order, then its errors in list order.
Eigen::VectorXcd sequential(const Circuit& circuit,
                            const std::vector<PauliError>& errors,
                            const Eigen::VectorXcd& initial) {
  StateVector state(initial);
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
  return state.to_vector();
}

Eigen::VectorXcd fused(const std::vector<FusedGate>& gates,
                       const Eigen::VectorXcd& initial) {
  StateVector state(initial);
  for (const FusedGate& gate : gates) {
    apply_fused(state, gate);
  }
  return state.to_vector();
}

Eigen::VectorXcd basis(int n, std::uint64_t index) {
  return Eigen::VectorXcd::Unit(Eigen::Index{1} << n,
                                static_cast<Eigen::Index>(index));
}

Eigen::VectorXcd random_state(int n, std::mt19937_64& rng) {
  std::normal_distribution<double> normal;
  Eigen::VectorXcd psi(Eigen::Index{1} << n);
  for (Eigen::Index i = 0; i < psi.size(); ++i) {
    psi(i) = cd(normal(rng), normal(rng));
  }
  return psi.normalized();
}

bool is_real(const FusedGate& gate) {
  return std::holds_alternative<Eigen::Matrix2d>(gate.matrix) ||
         std::holds_alternative<Eigen::Matrix4d>(gate.matrix);
}

bool is_one_qubit(const FusedGate& gate) {
  return std::holds_alternative<Eigen::Matrix2d>(gate.matrix) ||
         std::holds_alternative<Eigen::Matrix2cd>(gate.matrix);
}

// Depth-layer brickwork as brickwork_circuit builds it: RY on every qubit,
// then CNOTs on (0, 1), (2, 3), ... in odd layers and (1, 2), (3, 4), ... in
// even layers.
Circuit brickwork(int n, int depth, std::mt19937_64& rng) {
  std::uniform_real_distribution<double> angle(0, 2 * std::numbers::pi);
  Circuit circuit{n, {}};
  for (int layer = 1; layer <= depth; ++layer) {
    TimeStep rotations;
    for (int q = 0; q < n; ++q) {
      rotations.ops.push_back(Op{GateType::RY, {q, 0}, angle(rng)});
    }
    TimeStep entanglers;
    for (int q = layer % 2 == 1 ? 0 : 1; q + 1 < n; q += 2) {
      entanglers.ops.push_back(Op{GateType::CNOT, {q, q + 1}, 0.0});
    }
    circuit.timesteps.push_back(rotations);
    circuit.timesteps.push_back(entanglers);
  }
  return circuit;
}

// n_steps steps of up to 4 random gates each, with every gate type, random
// operands (both orders for two-qubit gates) and random angles. Qubits may
// repeat within a step: the IR allows it, and fusion must keep the order.
Circuit random_circuit(int n, int n_steps, std::mt19937_64& rng) {
  std::uniform_int_distribution<int> gate_index(0, kAllGates.size() - 1);
  std::uniform_int_distribution<int> qubit(0, n - 1);
  std::uniform_int_distribution<int> n_ops(0, 4);
  std::uniform_real_distribution<double> angle(0, 2 * std::numbers::pi);
  Circuit circuit{n, {}};
  for (int t = 0; t < n_steps; ++t) {
    TimeStep step;
    for (int i = n_ops(rng); i > 0; --i) {
      const GateType gate = kAllGates[gate_index(rng)];
      if (!is_two_qubit(gate)) {
        step.ops.push_back(Op{gate, {qubit(rng), 0}, angle(rng)});
      } else if (n >= 2) {
        const int a = qubit(rng);
        int b = qubit(rng);
        while (b == a) {
          b = qubit(rng);
        }
        step.ops.push_back(Op{gate, {a, b}, 0.0});
      }
    }
    circuit.timesteps.push_back(step);
  }
  return circuit;
}

// Up to 2 errors per step, sorted by timestep, with X, Y and Z.
std::vector<PauliError> random_errors(const Circuit& circuit,
                                      std::mt19937_64& rng) {
  constexpr std::array<GateType, 3> kPaulis = {GateType::X, GateType::Y,
                                               GateType::Z};
  std::uniform_int_distribution<int> qubit(0, circuit.n_qubits - 1);
  std::uniform_int_distribution<int> pauli(0, 2);
  std::uniform_int_distribution<int> count(0, 2);
  std::vector<PauliError> errors;
  for (std::size_t t = 0; t < circuit.timesteps.size(); ++t) {
    for (int i = count(rng); i > 0; --i) {
      errors.push_back(PauliError{t, qubit(rng), kPaulis[pauli(rng)]});
    }
  }
  return errors;
}

}  // namespace

TEST_CASE("fuse_trajectory applies each step's errors after its gates") {
  SECTION("H, then Z as an error of step 0, then H gives X") {
    // H Z H |0> = |1>. Z before the first H, or after the second, gives |0>.
    const TimeStep h0{{Op{GateType::H, {0, 0}, 0.0}}};
    const Circuit circuit{1, {h0, h0}};
    const std::vector<PauliError> errors{PauliError{0, 0, GateType::Z}};
    require_close(fused(fuse_trajectory(circuit, errors), basis(1, 0)),
                  basis(1, 1));
  }
  SECTION("CNOT 0 1, then X on qubit 0 as an error of the same step") {
    // CNOT then X0 gives |01>; X0 then CNOT would give |11>.
    const Circuit circuit{2, {TimeStep{{Op{GateType::CNOT, {0, 1}, 0.0}}}}};
    const std::vector<PauliError> errors{PauliError{0, 0, GateType::X}};
    require_close(fused(fuse_trajectory(circuit, errors), basis(2, 0)),
                  basis(2, 0b01));
  }
  SECTION("gates of a step keep their order") {
    // X0, then CNOT 0 1, in one step: |11>. The reverse order gives |01>.
    const Circuit circuit{2,
                          {TimeStep{{Op{GateType::X, {0, 0}, 0.0},
                                     Op{GateType::CNOT, {0, 1}, 0.0}}}}};
    require_close(fused(fuse_trajectory(circuit, {}), basis(2, 0)),
                  basis(2, 0b11));
  }
}

TEST_CASE("fused execution matches the general kernels on random circuits") {
  std::mt19937_64 rng(kSeed);
  for (const int n : {1, 2, 3, 5, 6}) {
    for (int trial = 0; trial < 20; ++trial) {
      CAPTURE(n, trial);
      const Circuit circuit = random_circuit(n, 12, rng);
      const std::vector<PauliError> errors = random_errors(circuit, rng);
      const Eigen::VectorXcd psi = random_state(n, rng);
      const std::vector<FusedGate> gates = fuse_trajectory(circuit, errors);
      for (const FusedGate& gate : gates) {
        if (!is_one_qubit(gate)) {
          REQUIRE(gate.qubits[0] < gate.qubits[1]);
        }
      }
      require_close(fused(gates, psi), sequential(circuit, errors, psi));
    }
  }
}

TEST_CASE("fuse_trajectory handles CNOT in both qubit orders") {
  std::mt19937_64 rng(kSeed);
  // Control below and above the target, adjacent and distant.
  const auto [control, target] = GENERATE(table<int, int>(
      {{0, 1}, {1, 0}, {0, 2}, {2, 0}, {1, 3}, {3, 1}, {0, 3}, {3, 0}}));
  CAPTURE(control, target);
  constexpr int kQubits = 4;
  // RY on both qubits before and after, so the CNOT fuses with 1q gates.
  const TimeStep rotations{{Op{GateType::RY, {control, 0}, 0.4},
                            Op{GateType::RY, {target, 0}, 1.3}}};
  const Circuit circuit{
      kQubits,
      {rotations, TimeStep{{Op{GateType::CNOT, {control, target}, 0.0}}},
       rotations}};
  const std::vector<FusedGate> gates = fuse_trajectory(circuit, {});
  REQUIRE(gates.size() == 1);
  CHECK(gates[0].qubits[0] == std::min(control, target));
  CHECK(gates[0].qubits[1] == std::max(control, target));
  const Eigen::VectorXcd psi = random_state(kQubits, rng);
  require_close(fused(gates, psi), sequential(circuit, {}, psi));
}

TEST_CASE("fuse_trajectory keeps real gates real") {
  const GateType gate = GENERATE(GateType::X, GateType::Z, GateType::H,
                                 GateType::RY, GateType::CNOT, GateType::CZ);
  CAPTURE(static_cast<int>(gate));
  const Op op =
      is_two_qubit(gate) ? Op{gate, {1, 0}, 0.0} : Op{gate, {0, 0}, 0.7};
  const Circuit circuit{2, {TimeStep{{op}}}};
  const std::vector<FusedGate> gates = fuse_trajectory(circuit, {});
  REQUIRE(gates.size() == 1);
  CHECK(is_real(gates[0]));
  CHECK(is_one_qubit(gates[0]) == !is_two_qubit(gate));
}

TEST_CASE("fuse_trajectory makes a cluster complex for a complex gate") {
  const GateType gate = GENERATE(GateType::Y, GateType::S, GateType::T,
                                 GateType::RX, GateType::RZ, GateType::CY);
  CAPTURE(static_cast<int>(gate));
  // Real H and CNOT, then the complex gate on the same two qubits.
  const Op op =
      is_two_qubit(gate) ? Op{gate, {1, 0}, 0.0} : Op{gate, {0, 0}, 0.7};
  const Circuit circuit{2,
                        {TimeStep{{Op{GateType::H, {0, 0}, 0.0},
                                   Op{GateType::CNOT, {0, 1}, 0.0}}},
                         TimeStep{{op}}}};
  const std::vector<FusedGate> gates = fuse_trajectory(circuit, {});
  REQUIRE(gates.size() == 1);
  CHECK_FALSE(is_real(gates[0]));
  CHECK_FALSE(is_one_qubit(gates[0]));
}

TEST_CASE(
    "fuse_trajectory: a Y error makes its cluster complex, X and Z do not") {
  const Circuit circuit{1, {TimeStep{{Op{GateType::H, {0, 0}, 0.0}}}}};
  for (const GateType pauli : {GateType::X, GateType::Y, GateType::Z}) {
    CAPTURE(static_cast<int>(pauli));
    const std::vector<FusedGate> gates =
        fuse_trajectory(circuit, {PauliError{0, 0, pauli}});
    REQUIRE(gates.size() == 1);
    CHECK(is_real(gates[0]) == (pauli != GateType::Y));
  }
}

TEST_CASE("Noiseless depth-20 brickwork fuses into 10 n - 1 real clusters") {
  std::mt19937_64 rng(kSeed);
  const int n = GENERATE(5, 12, 24);
  CAPTURE(n);
  const std::vector<FusedGate> gates =
      fuse_trajectory(brickwork(n, 20, rng), {});
  REQUIRE(gates.size() == static_cast<std::size_t>(10 * n - 1));
  std::size_t one_qubit = 0;
  for (const FusedGate& gate : gates) {
    CHECK(is_real(gate));
    one_qubit += is_one_qubit(gate) ? 1 : 0;
  }
  CHECK(one_qubit == 9);
}

TEST_CASE("Fused brickwork matches the general kernels") {
  std::mt19937_64 rng(kSeed);
  constexpr int kQubits = 7;
  const Circuit circuit = brickwork(kQubits, 20, rng);
  const Eigen::VectorXcd psi = random_state(kQubits, rng);
  require_close(fused(fuse_trajectory(circuit, {}), psi),
                sequential(circuit, {}, psi));
}

TEST_CASE("fuse_trajectory rejects bad errors with simulate's messages") {
  const Circuit circuit{3, {TimeStep{}, TimeStep{}}};
  CHECK_THROWS_AS(fuse_trajectory(circuit, {PauliError{2, 0, GateType::X}}),
                  std::out_of_range);
  CHECK_THROWS_WITH(fuse_trajectory(circuit, {PauliError{2, 0, GateType::X}}),
                    Equals("error timestep 2 out of range for 2 time steps"));
  CHECK_THROWS_AS(fuse_trajectory(circuit, {PauliError{1, 0, GateType::X},
                                            PauliError{0, 0, GateType::X}}),
                  std::invalid_argument);
  CHECK_THROWS_WITH(fuse_trajectory(circuit, {PauliError{1, 0, GateType::X},
                                              PauliError{0, 0, GateType::X}}),
                    Equals("errors must be sorted by timestep"));
  CHECK_THROWS_AS(fuse_trajectory(circuit, {PauliError{0, 0, GateType::H}}),
                  std::invalid_argument);
  CHECK_THROWS_WITH(fuse_trajectory(circuit, {PauliError{0, 0, GateType::H}}),
                    Equals("error pauli must be X, Y or Z"));
  CHECK_THROWS_AS(fuse_trajectory(circuit, {PauliError{0, 3, GateType::X}}),
                  std::out_of_range);
  CHECK_THROWS_WITH(fuse_trajectory(circuit, {PauliError{0, 3, GateType::X}}),
                    Equals("apply_1q: qubit 3 out of range for 3 qubits"));
}

TEST_CASE("fuse_trajectory rejects bad qubits with the kernels' messages") {
  // The messages that the per-gate kernels gave before fusion.
  const auto one_op = [](Op op) { return Circuit{3, {TimeStep{{op}}}}; };
  CHECK_THROWS_AS(fuse_trajectory(one_op(Op{GateType::RY, {3, 0}, 0.1}), {}),
                  std::out_of_range);
  CHECK_THROWS_WITH(fuse_trajectory(one_op(Op{GateType::RY, {3, 0}, 0.1}), {}),
                    Equals("apply_ry: qubit 3 out of range for 3 qubits"));
  CHECK_THROWS_WITH(fuse_trajectory(one_op(Op{GateType::H, {-1, 0}, 0.0}), {}),
                    Equals("apply_1q: qubit -1 out of range for 3 qubits"));
  CHECK_THROWS_WITH(
      fuse_trajectory(one_op(Op{GateType::CNOT, {0, 5}, 0.0}), {}),
      Equals("apply_cnot: qubit 5 out of range for 3 qubits"));
  CHECK_THROWS_AS(fuse_trajectory(one_op(Op{GateType::CNOT, {1, 1}, 0.0}), {}),
                  std::invalid_argument);
  CHECK_THROWS_WITH(
      fuse_trajectory(one_op(Op{GateType::CNOT, {1, 1}, 0.0}), {}),
      Equals("apply_cnot: both qubits are 1"));
  CHECK_THROWS_WITH(fuse_trajectory(one_op(Op{GateType::CZ, {4, 0}, 0.0}), {}),
                    Equals("apply_2q: qubit 4 out of range for 3 qubits"));
  CHECK_THROWS_WITH(fuse_trajectory(one_op(Op{GateType::CY, {2, 2}, 0.0}), {}),
                    Equals("apply_2q: both qubits are 2"));
  // The first operand is checked first.
  CHECK_THROWS_WITH(
      fuse_trajectory(one_op(Op{GateType::CNOT, {7, 5}, 0.0}), {}),
      Equals("apply_cnot: qubit 7 out of range for 3 qubits"));
  // qubits[1] of a one-qubit gate is unused, so it is not checked.
  CHECK_NOTHROW(fuse_trajectory(one_op(Op{GateType::X, {0, 9}, 0.0}), {}));
}

TEST_CASE("fuse_trajectory reports the first bad qubit in execution order") {
  // Step 0 has a bad error qubit, step 1 a bad gate qubit: the error runs
  // first.
  const Circuit circuit{
      2, {TimeStep{}, TimeStep{{Op{GateType::RY, {2, 0}, 0.1}}}}};
  CHECK_THROWS_WITH(fuse_trajectory(circuit, {PauliError{0, 4, GateType::Z}}),
                    Equals("apply_1q: qubit 4 out of range for 2 qubits"));
  // A gate runs before the errors of its own step.
  const Circuit same_step{2, {TimeStep{{Op{GateType::RY, {2, 0}, 0.1}}}}};
  CHECK_THROWS_WITH(fuse_trajectory(same_step, {PauliError{0, 4, GateType::Z}}),
                    Equals("apply_ry: qubit 2 out of range for 2 qubits"));
  // Error metadata is checked before any qubit.
  CHECK_THROWS_WITH(fuse_trajectory(same_step, {PauliError{0, 0, GateType::H}}),
                    Equals("error pauli must be X, Y or Z"));
}

TEST_CASE("apply_fused uses the kernel that matches the matrix") {
  std::mt19937_64 rng(kSeed);
  constexpr int kQubits = 3;
  const Eigen::VectorXcd psi = random_state(kQubits, rng);
  Eigen::Matrix4cd m4;
  Eigen::Matrix2cd m2;
  m4.setRandom();
  m2.setRandom();
  SECTION("complex 4x4: qubits[0] is the high bit") {
    StateVector expected(psi);
    expected.apply_2q(m4, 0, 2);
    StateVector state(psi);
    apply_fused(state, FusedGate{{0, 2}, m4});
    require_close(state.to_vector(), expected.to_vector());
  }
  SECTION("real 4x4") {
    StateVector expected(psi);
    expected.apply_2q(m4.real().cast<cd>(), 1, 2);
    StateVector state(psi);
    apply_fused(state, FusedGate{{1, 2}, Eigen::Matrix4d(m4.real())});
    require_close(state.to_vector(), expected.to_vector());
  }
  SECTION("complex 2x2") {
    StateVector expected(psi);
    expected.apply_1q(m2, 1);
    StateVector state(psi);
    apply_fused(state, FusedGate{{1, 0}, m2});
    require_close(state.to_vector(), expected.to_vector());
  }
  SECTION("real 2x2") {
    StateVector expected(psi);
    expected.apply_1q(m2.real().cast<cd>(), 2);
    StateVector state(psi);
    apply_fused(state, FusedGate{{2, 0}, Eigen::Matrix2d(m2.real())});
    require_close(state.to_vector(), expected.to_vector());
  }
}
