#ifndef ECM_QSIM_SIMULATOR_HPP_
#define ECM_QSIM_SIMULATOR_HPP_

// simulate(): one trajectory with given errors. run(): shots -> counts.

#include <vector>

#include "ecm_qsim/circuit.hpp"
#include "ecm_qsim/noise.hpp"
#include "ecm_qsim/state_vector.hpp"

namespace ecm_qsim {

// The state after one trajectory: starts from |0...0>, and for each time
// step t applies the gates of step t, then the errors with timestep t in
// list order. errors must be sorted by timestep, as sample_pauli_errors
// returns them.
// Throws std::out_of_range if an error's timestep is not a step of the
// circuit, std::invalid_argument if errors are not sorted by timestep or
// an error's pauli is not X, Y or Z. These checks run before allocation.
// StateVector errors propagate: std::invalid_argument for a bad
// n_qubits, std::length_error above StateVector::max_qubits(),
// std::out_of_range for a qubit outside the register.
StateVector simulate(const Circuit& circuit,
                     const std::vector<PauliError>& errors);

}  // namespace ecm_qsim

#endif  // ECM_QSIM_SIMULATOR_HPP_
