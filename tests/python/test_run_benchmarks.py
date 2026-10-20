import importlib.util
import json
import signal
import sys
from pathlib import Path

import pytest

from ecm_qsim import brickwork

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


MACHINE = {
    "cpu_model": "Test CPU", "physical_cores": 4, "memory_bytes": 8 << 30,
    "os": "TestOS-1", "gpu": "not used",
}  # fmt: skip


def run_args(tmp_path, label="reference", max_n=6):
    sim = tmp_path / "sim"
    sim.touch()
    return [
        "run", label, "--sim", str(sim), "--commit", "abc1234", "--max-n", str(max_n),
        "--peak-bandwidth-gbps", "100", "--compiler", "Apple clang 21",
        "--compiler-flags", "-O3 -DNDEBUG", "--results-dir", str(tmp_path / "results"),
    ]  # fmt: skip


def fake_measure(calls, capsys, max_n, fail_at=None):
    def measure(sim, circuit, *, n, shots, seed, noise):
        printed = capsys.readouterr().out
        calls.append({
            "n": n, "shots": shots, "seed": seed, "noise": noise,
            "circuit": Path(circuit).read_text(), "printed_before": printed,
        })  # fmt: skip
        if n == fail_at:
            raise RuntimeError(f"sim failed at n={n}")
        return {"counts": {"0" * n: shots}, "wall_s": 1e-6 * 2**n, "peak_rss_bytes": 16 << n}

    return measure


def test_run_measures_both_regimes_and_saves_result(tmp_path, monkeypatch, capsys):
    calls = []
    monkeypatch.setattr(rb, "machine_info", lambda: MACHINE)
    monkeypatch.setattr(rb, "measure", fake_measure(calls, capsys, max_n=6))
    assert rb.main(run_args(tmp_path)) == 0

    assert [c["n"] for c in calls] == [4, 5, 6, 12]
    assert all(c["seed"] == rb.SIM_SEED for c in calls)
    for c in calls[:3]:
        assert (c["shots"], c["noise"]) == (1, rb.NOISELESS)
        assert c["circuit"] == brickwork.brickwork_circuit(c["n"], 20, 42)
    noisy = calls[3]
    assert (noisy["shots"], noisy["noise"]) == (100_000, rb.MANY_SHOTS_NOISE)
    assert noisy["circuit"] == brickwork.brickwork_circuit(12, 20, 42)
    predicted = rb.predict_seconds(5, 1e-6 * 2**5, 6)
    assert "n=6" in calls[2]["printed_before"]
    assert f"{predicted:.3g}" in calls[2]["printed_before"]

    path = tmp_path / "results" / "reference.json"
    result = json.loads(path.read_text())
    assert result["label"] == "reference"
    assert result["commit"] == "abc1234"
    assert result["compiler"] == "Apple clang 21"
    assert result["compiler_flags"] == "-O3 -DNDEBUG"
    assert result["machine"] == MACHINE
    assert (result["depth"], result["angle_seed"], result["sim_seed"]) == (20, 42, 42)
    assert result["date"]

    points = result["many_qubits"]["points"]
    assert [p["n"] for p in points] == [4, 5, 6]
    bytes_6, flops_6 = rb.estimated_work(6)
    assert points[2]["modeled_bytes"] == bytes_6
    assert points[2]["useful_flops"] == flops_6
    assert points[2]["effective_bandwidth_gbps"] == pytest.approx(bytes_6 / 64e-6 / 1e9)
    assert points[2]["fraction_of_peak"] == pytest.approx(bytes_6 / 64e-6 / 1e9 / 100)
    assert result["many_qubits"]["prediction"] == {
        "n": 6, "from_n": 5, "predicted_s": pytest.approx(predicted),
        "measured_s": pytest.approx(64e-6),
    }  # fmt: skip

    shots = result["many_shots"]
    assert (shots["n"], shots["shots"]) == (12, 100_000)
    assert shots["noise"] == rb.MANY_SHOTS_NOISE
    assert shots["shots_per_s"] == pytest.approx(100_000 / (1e-6 * 2**12))
    assert str(path) in capsys.readouterr().out


def test_run_writes_nothing_when_a_measurement_fails(tmp_path, monkeypatch, capsys):
    monkeypatch.setattr(rb, "machine_info", lambda: MACHINE)
    monkeypatch.setattr(rb, "measure", fake_measure([], capsys, max_n=6, fail_at=12))
    with pytest.raises(RuntimeError, match="n=12"):
        rb.main(run_args(tmp_path))
    assert list((tmp_path / "results").glob("*")) == []


def test_run_rejects_max_n_without_a_point_to_predict_from(tmp_path, capsys):
    with pytest.raises(SystemExit) as exit_info:
        rb.main(run_args(tmp_path, max_n=4))
    assert exit_info.value.code == 2
    assert "--max-n must be at least 5" in capsys.readouterr().err
