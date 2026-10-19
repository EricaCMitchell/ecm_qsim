import subprocess
import sys
from pathlib import Path

import pytest

from ecm_qsim import sim

REPO_ROOT = Path(__file__).resolve().parents[2]


@pytest.mark.skipif(not sim.default_sim_path().exists(), reason="sim binary not built")
def test_script_writes_figure_and_prints_fit(tmp_path):
    out = tmp_path / "idle.png"
    result = subprocess.run(
        [
            sys.executable,
            "experiments/idle_noise.py",
            "--shots", "500",
            "--m", "0", "10", "40",
            "--p-idle", "0.02",
            "--out", str(out),
        ],  # fmt: skip
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
    assert out.stat().st_size > 0
    assert "p_idle fit" in result.stdout


@pytest.mark.skipif(not sim.default_sim_path().exists(), reason="sim binary not built")
def test_script_plots_p_idle_above_three_quarters_without_warnings(tmp_path):
    # Above 3/4 the decay base is negative, so the curves exist only at whole m.
    result = subprocess.run(
        [
            sys.executable,
            "-W", "error::RuntimeWarning",
            "experiments/idle_noise.py",
            "--shots", "500",
            "--m", "0", "1", "2", "3",
            "--p-idle", "0.9",
            "--out", str(tmp_path / "idle.png"),
        ],  # fmt: skip
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
