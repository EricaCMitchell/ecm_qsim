import importlib.util
import signal
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT = REPO_ROOT / "benchmarks" / "run_benchmarks.py"

_spec = importlib.util.spec_from_file_location("run_benchmarks", SCRIPT)
rb = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(rb)

MiB = 1 << 20


def test_run_child_reports_exit_code_stdout_and_peak_rss_in_bytes():
    child = "data = b'x' * (64 << 20); print(len(data))"
    code, out, err, wall_s, rss = rb.run_child([sys.executable, "-c", child])
    assert code == 0
    assert out == f"{64 * MiB}\n"
    assert err == ""
    assert wall_s > 0
    assert 64 * MiB <= rss < 1024 * MiB


def test_run_child_does_not_hang_on_large_stderr():
    child = "import sys; sys.stderr.write('e' * (1 << 20)); sys.exit(3)"
    code, _, err, _, _ = rb.run_child([sys.executable, "-c", child])
    assert code == 3
    assert len(err) == 1 << 20


def test_run_child_reports_signal_as_negative_code():
    child = "import os, signal; os.kill(os.getpid(), signal.SIGTERM)"
    code, *_ = rb.run_child([sys.executable, "-c", child])
    assert code == -signal.SIGTERM


def fake_child(code, stdout="", stderr="", wall_s=0.5, rss=1234, calls=None):
    def run_child(cmd):
        if calls is not None:
            calls.append(cmd)
        return code, stdout, stderr, wall_s, rss

    return run_child


def test_measure_builds_noiseless_command_and_returns_measurement(monkeypatch):
    calls = []
    stdout = '{"counts": {"0000": 1}}'
    monkeypatch.setattr(rb, "run_child", fake_child(0, stdout, calls=calls))
    result = rb.measure(
        Path("/x/sim"), Path("/c/n4.txt"), n=4, shots=1, seed=42, noise=rb.NOISELESS
    )
    assert calls == [[
        "/x/sim", "/c/n4.txt", "--shots", "1", "--seed", "42",
        "--p-gate", "0.0", "--p-idle", "0.0", "--p-readout", "0.0", "--json",
    ]]  # fmt: skip
    assert result == {"counts": {"0000": 1}, "wall_s": 0.5, "peak_rss_bytes": 1234}


def test_measure_builds_many_shots_command(monkeypatch):
    calls = []
    monkeypatch.setattr(rb, "run_child", fake_child(0, '{"counts": {}}', calls=calls))
    rb.measure(
        Path("sim"), Path("c.txt"), n=12, shots=rb.MANY_SHOTS, seed=rb.SIM_SEED,
        noise=rb.MANY_SHOTS_NOISE,
    )  # fmt: skip
    assert calls == [[
        "sim", "c.txt", "--shots", "100000", "--seed", "42",
        "--p-gate", "0.001", "--p-idle", "0.001", "--p-readout", "0.01", "--json",
    ]]  # fmt: skip


def test_measure_reports_nonzero_exit_status_with_n_and_stderr(monkeypatch):
    monkeypatch.setattr(rb, "run_child", fake_child(1, stderr="sim: bad circuit\n"))
    with pytest.raises(RuntimeError, match=r"(?s)n=7.*exit status 1.*sim: bad circuit"):
        rb.measure(Path("sim"), Path("c"), n=7, shots=1, seed=42, noise=rb.NOISELESS)


def test_measure_reports_terminating_signal(monkeypatch):
    monkeypatch.setattr(rb, "run_child", fake_child(-signal.SIGTERM, stderr="partial"))
    with pytest.raises(RuntimeError, match=r"(?s)n=29.*signal SIGTERM.*partial"):
        rb.measure(Path("sim"), Path("c"), n=29, shots=1, seed=42, noise=rb.NOISELESS)


def test_constants_match_required_regimes():
    assert (rb.DEPTH, rb.ANGLE_SEED, rb.SIM_SEED) == (20, 42, 42)
    assert (rb.MANY_SHOTS_N, rb.MANY_SHOTS) == (12, 100_000)
    assert rb.MANY_SHOTS_NOISE == {"p_gate": 1e-3, "p_idle": 1e-3, "p_readout": 1e-2}


def test_estimated_work_at_12_qubits():
    # 240 RY and 110 CNOT, each one pass over 2**12 amplitudes of 16 bytes,
    # read and written.
    assert rb.estimated_work(12) == (350 * 32 * 4096, 240 * 6 * 4096)


def test_predict_seconds_scales_by_modeled_bytes():
    # (240 + 110) * 2**12 / ((220 + 100) * 2**11) = 2.1875
    assert rb.predict_seconds(11, 2.0, 12) == pytest.approx(2 * 2.1875)


def fake_output(responses):
    def output(*cmd):
        return responses[cmd]

    return output


def test_machine_info_on_macos(monkeypatch):
    monkeypatch.setattr(rb.platform, "system", lambda: "Darwin")
    monkeypatch.setattr(rb, "_output", fake_output({
        ("sysctl", "-n", "machdep.cpu.brand_string", "hw.physicalcpu", "hw.memsize"):
            "Apple M3 Pro\n12\n25769803776\n",
    }))  # fmt: skip
    info = rb.machine_info()
    assert info["cpu_model"] == "Apple M3 Pro"
    assert info["physical_cores"] == 12
    assert info["memory_bytes"] == 25769803776
    assert info["gpu"] == "not used"
    assert info["os"]


def test_machine_info_on_linux(monkeypatch, tmp_path):
    meminfo = tmp_path / "meminfo"
    meminfo.write_text("MemTotal:       16318412 kB\nMemFree:  1000 kB\n")
    monkeypatch.setattr(rb, "MEMINFO", meminfo)
    monkeypatch.setattr(rb.platform, "system", lambda: "Linux")
    monkeypatch.setattr(rb, "_output", fake_output({
        ("lscpu",): "Architecture:  x86_64\nModel name:    AMD EPYC 7763 64-Core Processor\n",
        ("lscpu", "-p=CORE,SOCKET"): "# CORE,SOCKET\n0,0\n1,0\n0,0\n1,0\n0,1\n",
    }))  # fmt: skip
    info = rb.machine_info()
    assert info["cpu_model"] == "AMD EPYC 7763 64-Core Processor"
    assert info["physical_cores"] == 3
    assert info["memory_bytes"] == 16318412 * 1024
    assert info["gpu"] == "not used"


def test_machine_info_rejects_unsupported_platform(monkeypatch):
    monkeypatch.setattr(rb.platform, "system", lambda: "Windows")
    with pytest.raises(RuntimeError, match="Windows"):
        rb.machine_info()
