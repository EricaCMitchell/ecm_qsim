#ifndef ECM_QSIM_GATES_HPP_
#define ECM_QSIM_GATES_HPP_

// 2x2 and 4x4 gate matrices (Eigen fixed-size) per GateType.
// Grouped by matrix, not one file per gate: the kernels that apply them
// specialise on matrix shape (diagonal, permutation, dense), not gate name.

#include <Eigen/Core>

#include "ecm_qsim/circuit.hpp"

namespace ecm_qsim {

// Matrix of a one-qubit gate, acting on (amplitude with qubit = 0, amplitude
// with qubit = 1). theta is in radians and is used only by RX, RY and RZ:
// R_P(theta) = exp(-i theta P / 2). Fixed gates ignore it.
// Throws std::invalid_argument if gate is not a one-qubit gate.
Eigen::Matrix2cd matrix_1q(GateType gate, double theta = 0.0);

}  // namespace ecm_qsim

#endif  // ECM_QSIM_GATES_HPP_
