#ifndef ECM_QSIM_CIRCUIT_HPP_
#define ECM_QSIM_CIRCUIT_HPP_

// Intermediate representation of a parsed circuit.
// Time steps are explicit because idle noise depends on which qubits are idle
// in each TICK.

#include <array>
#include <vector>

namespace ecm_qsim {

enum class GateType { X, Y, Z, H, S, T, RX, RY, RZ, CNOT, CZ, CY };

struct Op {
  GateType gate;
  // Control first for two-qubit gates; [1] unused for one-qubit gates.
  std::array<int, 2> qubits;
  double theta;  // rotation angle in radians; unused for fixed gates
};

struct TimeStep {
  std::vector<Op> ops;  // one TICK; may be empty
};

struct Circuit {
  int n_qubits;
  std::vector<TimeStep> timesteps;
};

}  // namespace ecm_qsim

#endif  // ECM_QSIM_CIRCUIT_HPP_
