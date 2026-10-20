#ifndef ECM_QSIM_FUSION_HPP_
#define ECM_QSIM_FUSION_HPP_

// Gate fusion: one trajectory (circuit gates and sampled Pauli errors) as
// fused gates on at most two qubits, so that each pass over the state
// applies several gates.

#include <Eigen/Core>
#include <array>
#include <variant>
#include <vector>

#include "ecm_qsim/circuit.hpp"
#include "ecm_qsim/noise.hpp"
#include "ecm_qsim/state_vector.hpp"

namespace ecm_qsim {

// The alternative gives the arity (2x2: one qubit, 4x4: two qubits) and
// whether the matrix is real.
using FusedMatrix = std::variant<Eigen::Matrix2d, Eigen::Matrix4d,
                                 Eigen::Matrix2cd, Eigen::Matrix4cd>;

struct FusedGate {
  // One-qubit matrix: qubits[0]; qubits[1] is 0 and unused.
  // Two-qubit matrix: qubits[0] < qubits[1], and qubits[0] is the high bit
  // of the matrix index, as j is in StateVector::apply_2q.
  std::array<int, 2> qubits;
  FusedMatrix matrix;
};

// The checks that simulate() runs on errors before it allocates the state.
// Throws std::out_of_range if an error's timestep is not a step of the
// circuit, std::invalid_argument if errors are not sorted by timestep or
// an error's pauli is not X, Y or Z.
void validate_pauli_errors(const Circuit& circuit,
                           const std::vector<PauliError>& errors);

// The trajectory of simulate(): for each step, its gates in order, then its
// errors in list order, as fused gates whose product is the same operator.
// Greedy fusion to width two: each operation joins the open clusters on its
// qubits if their union has at most two qubits. Otherwise those clusters
// are emitted and the operation starts a new cluster. Only clusters on
// disjoint qubits change order.
// Throws as validate_pauli_errors does, then, in execution order, the
// exception and message of the kernel that simulate() used before fusion
// for the first operation with a bad qubit (for example "apply_ry: qubit 3
// out of range for 3 qubits").
std::vector<FusedGate> fuse_trajectory(const Circuit& circuit,
                                       const std::vector<PauliError>& errors);

// Applies gate with apply_real_1q, apply_real_2q, apply_1q or apply_2q, as
// its matrix alternative selects.
void apply_fused(StateVector& state, const FusedGate& gate);

}  // namespace ecm_qsim

#endif  // ECM_QSIM_FUSION_HPP_
