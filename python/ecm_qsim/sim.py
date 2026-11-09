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

"""Runs the C++ `sim` binary as a subprocess."""

import json
import os
import subprocess
import tempfile
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parents[2]


class SimError(RuntimeError):
    """The sim binary is missing or exited with an error."""


def default_sim_path() -> Path:
    """$ECM_QSIM_SIM if set, else build/sim in the repository."""
    return Path(os.environ.get("ECM_QSIM_SIM", _REPO_ROOT / "build" / "sim"))


def run(
    circuit: str,
    *,
    shots: int,
    seed: int,
    p_gate: float = 0.0,
    p_idle: float = 0.0,
    p_readout: float = 0.0,
    sim_path: Path | None = None,
) -> dict[str, int]:
    """Runs the circuit text and returns the counts per measured bitstring."""
    sim_path = default_sim_path() if sim_path is None else Path(sim_path)
    if not sim_path.exists():
        raise SimError(
            f"sim binary not found at {sim_path}; build it with "
            "'cmake -S . -B build -G Ninja && cmake --build build'"
        )
    with tempfile.TemporaryDirectory() as tmp:
        circuit_path = Path(tmp) / "circuit.txt"
        circuit_path.write_text(circuit)
        # repr is the shortest form that round-trips, so no precision is lost.
        args = [
            str(sim_path),
            str(circuit_path),
            "--shots", str(shots),
            "--seed", str(seed),
            "--p-gate", repr(float(p_gate)),
            "--p-idle", repr(float(p_idle)),
            "--p-readout", repr(float(p_readout)),
            "--json",
        ]  # fmt: skip
        result = subprocess.run(args, capture_output=True, text=True)
    if result.returncode != 0:
        raise SimError(result.stderr.strip())
    return json.loads(result.stdout)["counts"]
