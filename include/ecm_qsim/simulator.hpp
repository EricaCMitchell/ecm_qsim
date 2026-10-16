#ifndef ECM_QSIM_SIMULATOR_HPP_
#define ECM_QSIM_SIMULATOR_HPP_

// simulate(): one trajectory with given errors. run(): shots -> counts.

#include <cstdint>
#include <map>
#include <string>
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

// Measurement histogram: n-character bitstring (most significant bit
// first, qubit 0 rightmost) -> number of shots. Sorted, so the order of
// outcomes is fixed.
using Counts = std::map<std::string, int>;

// Runs shots Monte Carlo trajectories and counts the measured
// bitstrings. One Rng is seeded with seed. Each shot draws from it in
// this order: sample_pauli_errors, one draw for the measurement
// (uniform_from_bits), flip_readout_bits. Every shot starts from
// |0...0>.
// Throws std::invalid_argument if shots < 1 or a noise probability is not
// in [0, 1], both checked before any work. Errors from sample_pauli_errors
// and simulate propagate (std::length_error if the register does not fit
// in memory).
Counts run(const Circuit& circuit, const NoiseModel& noise, int shots,
           std::uint64_t seed);

}  // namespace ecm_qsim

#endif  // ECM_QSIM_SIMULATOR_HPP_
