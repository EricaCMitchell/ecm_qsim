import re
import subprocess
import sys
from pathlib import Path

from ecm_qsim import brickwork

REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT = REPO_ROOT / "benchmarks" / "brickwork_generator.py"


def run_script(*args: str, cwd: Path = REPO_ROOT) -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, str(SCRIPT), *args],
        cwd=cwd,
        capture_output=True,
        text=True,
    )


def test_script_writes_circuit_to_out(tmp_path):
    out = tmp_path / "c.txt"
    result = run_script(
        "--qubits", "4", "--depth", "3", "--seed", "5", "--out", str(out)
    )  # fmt: skip
    assert result.returncode == 0, result.stderr
    assert out.read_text() == brickwork.brickwork_circuit(4, 3, 5)
    assert str(out) in result.stdout


def test_script_without_seed_gives_new_angles_each_run(tmp_path):
    out = tmp_path / "c.txt"
    texts = []
    for _ in range(2):
        result = run_script("--qubits", "4", "--depth", "3", "--out", str(out))
        assert result.returncode == 0, result.stderr
        text = out.read_text()
        header = re.fullmatch(r"# brickwork n=4 depth=3 seed=(\d+)", text.splitlines()[0])
        seed = int(header[1])
        assert text == brickwork.brickwork_circuit(4, 3, seed)
        texts.append(text)
    assert texts[0] != texts[1]


def test_script_default_path_ignores_working_directory(tmp_path):
    expected = REPO_ROOT / "benchmarks" / "circuits" / "brickwork_n3_d1.txt"
    try:
        result = run_script("--qubits", "3", "--depth", "1", "--seed", "1", cwd=tmp_path)
        assert result.returncode == 0, result.stderr
        assert expected.read_text() == brickwork.brickwork_circuit(3, 1, 1)
        assert list(tmp_path.iterdir()) == []
    finally:
        expected.unlink(missing_ok=True)


def test_script_rejects_invalid_qubits(tmp_path):
    out = tmp_path / "c.txt"
    result = run_script("--qubits", "1", "--out", str(out))
    assert result.returncode == 2
    assert "qubits" in result.stderr
    assert not out.exists()
