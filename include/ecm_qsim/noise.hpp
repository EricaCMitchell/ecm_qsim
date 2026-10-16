#ifndef ECM_QSIM_NOISE_HPP_
#define ECM_QSIM_NOISE_HPP_

// Gate, idle and readout noise: trajectory sampling and readout flips.

#include <cstddef>
#include <cstdint>
#include <random>
#include <vector>

#include "ecm_qsim/circuit.hpp"

namespace ecm_qsim {

struct NoiseModel {
  double p_gate = 0.0;
  double p_idle = 0.0;
  double p_readout = 0.0;
};

// Throws std::invalid_argument("<field> must be in [0, 1], got <value>") for
// the first field, in the order p_gate, p_idle, p_readout, that is not in
// [0, 1] (NaN included).
void validate_noise_model(const NoiseModel& noise);

// The caller-owned RNG for noise and measurement draws.
using Rng = std::mt19937_64;

// A Pauli error on one qubit, applied after all gates of time step
// `timestep` (an index into Circuit::timesteps).
struct PauliError {
  std::size_t timestep;
  int qubit;
  GateType pauli;  // X, Y or Z
  bool operator==(const PauliError&) const = default;
};

// Samples the gate and idle errors of one shot, sorted by (timestep,
// qubit). Throws std::invalid_argument if a probability is not in [0, 1]
// or circuit.n_qubits is not in 1..kMaxQubits. Throws std::out_of_range
// if an operation has a qubit outside the circuit. All validation checks
// run before the first draw, so validation failures leave rng unchanged.
std::vector<PauliError> sample_pauli_errors(const Circuit& circuit,
                                            const NoiseModel& noise, Rng& rng);

// Flips each of the n_qubits low bits of index independently with
// probability p_readout, drawing for qubit 0 first. No draws if p_readout
// is 0. Throws std::invalid_argument if p_readout is not in [0, 1] or
// n_qubits is not in 1..kMaxQubits, std::out_of_range if
// index >= 2^n_qubits.
std::uint64_t flip_readout_bits(std::uint64_t index, int n_qubits,
                                double p_readout, Rng& rng);

}  // namespace ecm_qsim

#endif  // ECM_QSIM_NOISE_HPP_
