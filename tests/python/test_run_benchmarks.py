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


def synthetic_result(label, walls, shots_per_s):
    largest = max(walls)
    result = {
        "label": label, "sim": f"/bin/{label}", "commit": "abc", "date": "2026-10-19",
        "machine": dict(MACHINE), "compiler": "Apple clang 21",
        "compiler_flags": "-O3 -DNDEBUG", "peak_bandwidth_gbps": 100.0,
        "depth": 20, "angle_seed": 42, "sim_seed": 42,
        "many_qubits": {
            "points": [
                {
                    "n": n, "wall_s": s, "peak_rss_bytes": 16 << n,
                    "modeled_bytes": rb.estimated_work(n)[0],
                    "useful_flops": rb.estimated_work(n)[1],
                    "effective_bandwidth_gbps": rb.estimated_work(n)[0] / s / 1e9,
                    "fraction_of_peak": rb.estimated_work(n)[0] / s / 1e9 / 100,
                }
                for n, s in sorted(walls.items())
            ],
            "prediction": {
                "n": largest, "from_n": largest - 1,
                "predicted_s": 0.9 * walls[largest], "measured_s": walls[largest],
            },
        },
        "many_shots": {
            "n": 12, "shots": 100_000, "noise": rb.MANY_SHOTS_NOISE,
            "wall_s": 100_000 / shots_per_s, "peak_rss_bytes": 1 << 20,
            "shots_per_s": shots_per_s,
        },
    }  # fmt: skip
    return result


REFERENCE = synthetic_result("reference", {4: 1.0, 5: 2.0, 6: 4.0}, 1000.0)
OPTIMIZED = synthetic_result("optimized", {4: 0.5, 5: 0.5, 6: 1.0, 7: 2.0}, 4000.0)


def write_results(directory, *results):
    directory.mkdir(parents=True, exist_ok=True)
    for result in results:
        (directory / f"{result['label']}.json").write_text(json.dumps(result))


def test_draw_plots_one_log_scaled_line_per_version():
    figure = rb.draw([REFERENCE, OPTIMIZED])
    (axes,) = figure.axes
    assert axes.get_yscale() == "log"
    assert axes.get_xlabel() == "qubits n"
    assert [line.get_label() for line in axes.lines if not line.get_label().startswith("_")][
        :2
    ] == ["reference", "optimized"]
    reference_line = next(line for line in axes.lines if line.get_label() == "reference")
    assert list(reference_line.get_xdata()) == [4, 5, 6]
    assert list(reference_line.get_ydata()) == [1.0, 2.0, 4.0]


def test_draw_marks_each_prediction_with_a_hollow_marker():
    figure = rb.draw([REFERENCE, OPTIMIZED])
    predictions = [
        line for line in figure.axes[0].lines
        if line.get_marker() == "o" and line.get_markerfacecolor() == "none"
    ]  # fmt: skip
    assert sorted((list(p.get_xdata()), list(p.get_ydata())) for p in predictions) == [
        ([6], [pytest.approx(3.6)]),
        ([7], [pytest.approx(1.8)]),
    ]


def test_summary_reports_speed_ups_against_the_reference():
    reference, optimized = rb.summary([REFERENCE, OPTIMIZED])
    assert reference["label"] == "reference"
    assert (reference["largest_n"], reference["many_qubits_speed_up"]) == (6, 1.0)
    assert optimized["largest_n"] == 7
    assert optimized["peak_rss_bytes"] == 16 << 7
    assert optimized["effective_bandwidth_gbps"] == pytest.approx(
        rb.estimated_work(7)[0] / 2.0 / 1e9
    )
    assert optimized["fraction_of_peak"] == pytest.approx(rb.estimated_work(7)[0] / 2e11)
    assert optimized["shots_per_s"] == 4000.0
    # Largest common n is 6: 4.0 s against 1.0 s.
    assert optimized["common_n"] == 6
    assert optimized["many_qubits_speed_up"] == pytest.approx(4.0)
    assert optimized["many_shots_speed_up"] == pytest.approx(4.0)


def test_plot_writes_png_and_prints_table(tmp_path, capsys):
    write_results(tmp_path, REFERENCE, OPTIMIZED)
    assert rb.main(["plot", "--results-dir", str(tmp_path)]) == 0
    png = tmp_path / "time_vs_n.png"
    assert png.stat().st_size > 0
    captured = capsys.readouterr()
    assert str(png) in captured.out
    assert "optimized" in captured.out and "4.00x" in captured.out


def test_plot_requires_a_reference_result(tmp_path):
    write_results(tmp_path, OPTIMIZED)
    with pytest.raises(SystemExit, match="reference.json"):
        rb.main(["plot", "--results-dir", str(tmp_path)])
