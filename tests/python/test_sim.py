# Copyright 2026 Erica C. Mitchell
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
# SPDX-License-Identifier: Apache-2.0

import pytest

from ecm_qsim import sim

needs_sim = pytest.mark.skipif(
    not sim.default_sim_path().exists(), reason="sim binary not built"
)


@needs_sim
def test_run_returns_counts_for_circuit_text():
    assert sim.run("QUBITS 2\nX 1\n", shots=50, seed=1) == {"10": 50}


@needs_sim
def test_run_passes_noise_parameters():
    assert sim.run("QUBITS 1\nX 0\n", shots=100, seed=1, p_readout=1.0) == {"0": 100}


@needs_sim
def test_run_is_reproducible_for_a_seed():
    circuit = "QUBITS 1\nH 0\n"
    assert sim.run(circuit, shots=200, seed=3) == sim.run(circuit, shots=200, seed=3)


@needs_sim
def test_run_reports_sim_errors():
    with pytest.raises(sim.SimError, match="unknown gate"):
        sim.run("QUBITS 1\nFOO 0\n", shots=1, seed=1)


def test_run_with_missing_binary_says_how_to_build(tmp_path):
    with pytest.raises(sim.SimError, match="cmake --build build"):
        sim.run("QUBITS 1\n", shots=1, seed=1, sim_path=tmp_path / "nope")
