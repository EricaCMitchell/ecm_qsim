#ifndef ECM_QSIM_CIRCUIT_HPP_
#define ECM_QSIM_CIRCUIT_HPP_

// Intermediate representation of a parsed circuit.
// Time steps are explicit because idle noise depends on which qubits are idle
// in each TICK.

#include <array>
#include <vector>

namespace ecm_qsim {

enum class GateType { X, Y, Z, H, S, T, RX, RY, RZ, CNOT, CZ, CY };

// True for the controlled gates, which act on qubits[0] and qubits[1]. No
// default case, so -Wswitch flags a gate added later.
constexpr bool is_two_qubit(GateType gate) noexcept {
  switch (gate) {
    case GateType::X:
    case GateType::Y:
    case GateType::Z:
    case GateType::H:
    case GateType::S:
    case GateType::T:
    case GateType::RX:
    case GateType::RY:
    case GateType::RZ:
      return false;
    case GateType::CNOT:
    case GateType::CZ:
    case GateType::CY:
      return true;
  }
  return false;
}

struct Op {
  GateType gate;
  // Control first for two-qubit gates; [1] is 0 and unused for one-qubit
  // gates.
  std::array<int, 2> qubits;
  double theta;  // rotation angle in radians; 0 and unused for fixed gates

  bool operator==(const Op&) const = default;
};

struct TimeStep {
  std::vector<Op> ops;  // one TICK; may be empty
};

// Largest qubit count the IR represents: the amplitude count 2^n fits in
// uint64_t. Whether the state vector fits in memory is checked where it is
// allocated, not here.
inline constexpr int kMaxQubits = 63;

struct Circuit {
  int n_qubits;  // 1..kMaxQubits
  std::vector<TimeStep> timesteps;
};

}  // namespace ecm_qsim

#endif  // ECM_QSIM_CIRCUIT_HPP_
