// Copyright 2026 Erica C. Mitchell
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//
// SPDX-License-Identifier: Apache-2.0

#include "ecm_qsim/simulator.hpp"

#include <stdexcept>
#include <string>

#include "ecm_qsim/fusion.hpp"
#include "ecm_qsim/measurement.hpp"

namespace ecm_qsim {

StateVector simulate(const Circuit& circuit,
                     const std::vector<PauliError>& errors) {
  // Checked before allocation. fuse_trajectory repeats these checks.
  validate_pauli_errors(circuit, errors);
  StateVector state(circuit.n_qubits);
  for (const FusedGate& gate : fuse_trajectory(circuit, errors)) {
    apply_fused(state, gate);
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
    // One full expression, so the temporary state outlives the sampling.
    const std::uint64_t index =
        sample_index(simulate(circuit, errors), uniform_from_bits(rng()));
    const std::uint64_t measured =
        flip_readout_bits(index, circuit.n_qubits, noise.p_readout, rng);
    ++counts[to_bitstring(measured, circuit.n_qubits)];
  }
  return counts;
}

}  // namespace ecm_qsim
