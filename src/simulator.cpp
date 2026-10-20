#include "ecm_qsim/simulator.hpp"

#include <cstddef>
#include <stdexcept>
#include <string>

#include "ecm_qsim/gates.hpp"
#include "ecm_qsim/measurement.hpp"

namespace ecm_qsim {

namespace {

// RY and CNOT, the brickwork benchmark's gates, take the specialized
// kernels; every other gate takes the general matrix kernels.
void apply_op(StateVector& state, const Op& op) {
  if (op.gate == GateType::RY) {
    state.apply_ry(op.theta, op.qubits[0]);
  } else if (op.gate == GateType::CNOT) {
    state.apply_cnot(op.qubits[0], op.qubits[1]);
  } else if (is_two_qubit(op.gate)) {
    state.apply_2q(matrix_2q(op.gate), op.qubits[0], op.qubits[1]);
  } else {
    state.apply_1q(matrix_1q(op.gate, op.theta), op.qubits[0]);
  }
}

}  // namespace

StateVector simulate(const Circuit& circuit,
                     const std::vector<PauliError>& errors) {
  const std::size_t n_steps = circuit.timesteps.size();
  for (std::size_t i = 0; i < errors.size(); ++i) {
    const PauliError& error = errors[i];
    if (error.timestep >= n_steps) {
      throw std::out_of_range(
          "error timestep " + std::to_string(error.timestep) +
          " out of range for " + std::to_string(n_steps) + " time steps");
    }
    if (i > 0 && error.timestep < errors[i - 1].timestep) {
      throw std::invalid_argument("errors must be sorted by timestep");
    }
    if (error.pauli != GateType::X && error.pauli != GateType::Y &&
        error.pauli != GateType::Z) {
      throw std::invalid_argument("error pauli must be X, Y or Z");
    }
  }

  StateVector state(circuit.n_qubits);
  std::size_t next = 0;
  for (std::size_t t = 0; t < n_steps; ++t) {
    for (const Op& op : circuit.timesteps[t].ops) {
      apply_op(state, op);
    }
    for (; next < errors.size() && errors[next].timestep == t; ++next) {
      state.apply_1q(matrix_1q(errors[next].pauli), errors[next].qubit);
    }
  }
  return state;
}

Counts run(const Circuit& circuit, const NoiseModel& noise, int shots,
           std::uint64_t seed) {
  if (shots < 1) {
    throw std::invalid_argument("shots must be at least 1, got " +
                                std::to_string(shots));
  }
  validate_noise_model(noise);

  Rng rng(seed);
  Counts counts;
  for (int shot = 0; shot < shots; ++shot) {
    const std::vector<PauliError> errors =
        sample_pauli_errors(circuit, noise, rng);
    const std::uint64_t index = sample_index(
        simulate(circuit, errors).probabilities(), uniform_from_bits(rng()));
    const std::uint64_t measured =
        flip_readout_bits(index, circuit.n_qubits, noise.p_readout, rng);
    ++counts[to_bitstring(measured, circuit.n_qubits)];
  }
  return counts;
}

}  // namespace ecm_qsim
