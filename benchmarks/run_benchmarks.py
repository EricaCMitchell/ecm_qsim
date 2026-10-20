"""Time a sim binary on the brickwork benchmark circuits and plot the results.

`run` measures one named version in the two required regimes and writes
benchmarks/results/<LABEL>.json. `plot` reads every result in the directory
and draws wall time against qubit count.
"""

import argparse
import json
import os
import platform
import signal
import subprocess
import sys
import tempfile
import time
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.ticker import MaxNLocator  # noqa: E402

from ecm_qsim import brickwork  # noqa: E402

DEPTH = 20
ANGLE_SEED = 42
SIM_SEED = 42
MANY_SHOTS_N = 12
MANY_SHOTS = 100_000
NOISELESS = {"p_gate": 0.0, "p_idle": 0.0, "p_readout": 0.0}
MANY_SHOTS_NOISE = {"p_gate": 1e-3, "p_idle": 1e-3, "p_readout": 1e-2}

# One gate reads and writes every complex<double> amplitude once.
BYTES_PER_AMPLITUDE = 32
RY_FLOPS_PER_AMPLITUDE = 6

# Categorical colours in fixed order; the reference always takes the first.
COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300"]

RESULTS_DIR = Path(__file__).resolve().parent / "results"
MEMINFO = Path("/proc/meminfo")


def run_child(cmd: list[str]) -> tuple[int, str, str, float, int]:
    """Runs cmd and returns exit code, stdout, stderr, wall seconds, peak RSS bytes.

    The exit code is negative when a signal killed the child. Output goes to
    temporary files so that a full pipe cannot block the child.
    """
    with tempfile.TemporaryFile() as out, tempfile.TemporaryFile() as err:
        start = time.perf_counter()
        proc = subprocess.Popen(cmd, stdout=out, stderr=err)
        # wait4 reaps the child itself; Popen.wait() would discard its usage.
        _, status, usage = os.wait4(proc.pid, 0)
        wall_s = time.perf_counter() - start
        proc.returncode = os.waitstatus_to_exitcode(status)
        out.seek(0)
        err.seek(0)
        stdout = out.read().decode()
        stderr = err.read().decode()
    # ru_maxrss is in bytes on macOS and in kilobytes on Linux.
    rss_unit = 1 if platform.system() == "Darwin" else 1024
    return proc.returncode, stdout, stderr, wall_s, usage.ru_maxrss * rss_unit


def estimated_work(n: int, depth: int = DEPTH) -> tuple[int, int]:
    """Modeled bytes moved and useful FLOPs for the n-qubit brickwork circuit."""
    ry = depth * n
    odd_layers = (depth + 1) // 2
    cnot = odd_layers * (n // 2) + (depth - odd_layers) * ((n - 1) // 2)
    amplitudes = 2**n
    bytes_moved = (ry + cnot) * BYTES_PER_AMPLITUDE * amplitudes
    flops = ry * RY_FLOPS_PER_AMPLITUDE * amplitudes
    return bytes_moved, flops


def predict_seconds(previous_n: int, previous_s: float, n: int) -> float:
    """Scales the time measured at previous_n by the ratio of modeled bytes."""
    return previous_s * estimated_work(n)[0] / estimated_work(previous_n)[0]


def _output(*cmd: str) -> str:
    return subprocess.run(cmd, capture_output=True, text=True, check=True).stdout


def machine_info() -> dict[str, object]:
    """CPU model, physical cores, memory and OS of this machine."""
    system = platform.system()
    if system == "Darwin":
        cpu_model, cores, memory = _output(
            "sysctl", "-n", "machdep.cpu.brand_string", "hw.physicalcpu", "hw.memsize"
        ).splitlines()
        physical_cores = int(cores)
        memory_bytes = int(memory)
    elif system == "Linux":
        cpu_model = next(
            line.split(":", 1)[1]
            for line in _output("lscpu").splitlines()
            if line.startswith("Model name:")
        )
        physical_cores = len({
            line
            for line in _output("lscpu", "-p=CORE,SOCKET").splitlines()
            if line and not line.startswith("#")
        })  # fmt: skip
        mem_total = next(
            line for line in MEMINFO.read_text().splitlines() if line.startswith("MemTotal:")
        )
        memory_bytes = int(mem_total.split()[1]) * 1024
    else:
        raise RuntimeError(f"unsupported platform {system!r}: only macOS and Linux")
    return {
        "cpu_model": cpu_model.strip(),
        "physical_cores": physical_cores,
        "memory_bytes": memory_bytes,
        "os": platform.platform(),
        "gpu": "not used",
    }


def measure(
    sim: Path, circuit: Path, *, n: int, shots: int, seed: int, noise: dict[str, float]
) -> dict[str, object]:
    """Runs sim once on circuit and returns its counts, wall time and peak RSS."""
    cmd = [
        str(sim), str(circuit),
        "--shots", str(shots),
        "--seed", str(seed),
        "--p-gate", repr(float(noise["p_gate"])),
        "--p-idle", repr(float(noise["p_idle"])),
        "--p-readout", repr(float(noise["p_readout"])),
        "--json",
    ]  # fmt: skip
    code, stdout, stderr, wall_s, rss = run_child(cmd)
    if code != 0:
        how = (
            f"killed by signal {signal.Signals(-code).name}"
            if code < 0
            else f"exit status {code}"
        )
        raise RuntimeError(f"sim failed at n={n}: {how}\n{stderr.strip()}")
    return {"counts": json.loads(stdout)["counts"], "wall_s": wall_s, "peak_rss_bytes": rss}


def run(args: argparse.Namespace) -> None:
    """Measures both regimes for one version and writes <label>.json."""
    result = {
        "label": args.label,
        "sim": str(args.sim.resolve()),
        "commit": args.commit,
        "date": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "machine": machine_info(),
        "compiler": args.compiler,
        "compiler_flags": args.compiler_flags,
        "peak_bandwidth_gbps": args.peak_bandwidth_gbps,
        "depth": DEPTH,
        "angle_seed": ANGLE_SEED,
        "sim_seed": SIM_SEED,
    }
    with tempfile.TemporaryDirectory() as tmp:

        def circuit(n: int) -> Path:
            path = Path(tmp) / f"brickwork_n{n}_d{DEPTH}.txt"
            path.write_text(brickwork.brickwork_circuit(n, DEPTH, ANGLE_SEED))
            return path

        points = []
        for n in range(4, args.max_n + 1):
            if n == args.max_n:
                previous = points[-1]
                predicted_s = predict_seconds(previous["n"], previous["wall_s"], n)
                print(
                    f"predicted n={n}: {predicted_s:.3g} s "
                    f"(from n={previous['n']}: {previous['wall_s']:.3g} s)",
                    flush=True,
                )
            m = measure(args.sim, circuit(n), n=n, shots=1, seed=SIM_SEED, noise=NOISELESS)
            modeled_bytes, useful_flops = estimated_work(n)
            bandwidth_gbps = modeled_bytes / m["wall_s"] / 1e9
            points.append({
                "n": n,
                "wall_s": m["wall_s"],
                "peak_rss_bytes": m["peak_rss_bytes"],
                "modeled_bytes": modeled_bytes,
                "useful_flops": useful_flops,
                "effective_bandwidth_gbps": bandwidth_gbps,
                "fraction_of_peak": bandwidth_gbps / args.peak_bandwidth_gbps,
            })  # fmt: skip
            print(f"n={n}: {m['wall_s']:.3g} s", flush=True)
        result["many_qubits"] = {
            "points": points,
            "prediction": {
                "n": args.max_n,
                "from_n": args.max_n - 1,
                "predicted_s": predicted_s,
                "measured_s": points[-1]["wall_s"],
            },
        }

        m = measure(
            args.sim, circuit(MANY_SHOTS_N), n=MANY_SHOTS_N, shots=MANY_SHOTS,
            seed=SIM_SEED, noise=MANY_SHOTS_NOISE,
        )  # fmt: skip
        result["many_shots"] = {
            "n": MANY_SHOTS_N,
            "shots": MANY_SHOTS,
            "noise": MANY_SHOTS_NOISE,
            "wall_s": m["wall_s"],
            "peak_rss_bytes": m["peak_rss_bytes"],
            "shots_per_s": MANY_SHOTS / m["wall_s"],
        }

    args.results_dir.mkdir(parents=True, exist_ok=True)
    path = args.results_dir / f"{args.label}.json"
    partial = path.with_name(f".{path.name}.tmp")
    partial.write_text(json.dumps(result, indent=2) + "\n")
    os.replace(partial, path)

    largest = points[-1]
    print(path)
    print(
        f"n={args.max_n}: predicted {predicted_s:.3g} s, measured {largest['wall_s']:.3g} s, "
        f"peak RSS {largest['peak_rss_bytes'] / 2**30:.3g} GiB"
    )
    print(f"n={MANY_SHOTS_N}, {MANY_SHOTS} noisy shots: {MANY_SHOTS / m['wall_s']:.4g} shots/s")


def load_results(results_dir: Path) -> list[dict]:
    """Every result in results_dir, with the reference first."""
    results = {
        path.stem: json.loads(path.read_text()) for path in sorted(results_dir.glob("*.json"))
    }
    if "reference" not in results:
        sys.exit(f"no reference.json in {results_dir}: run the reference first")
    reference = results.pop("reference")
    return [reference, *results.values()]


def draw(results: list[dict]) -> plt.Figure:
    """Wall time against n, one line per version, predictions as hollow markers."""
    figure, axes = plt.subplots(figsize=(7, 4.5), layout="constrained")
    for result, color in zip(results, COLORS):
        points = result["many_qubits"]["points"]
        axes.plot(
            [p["n"] for p in points], [p["wall_s"] for p in points],
            color=color, linewidth=2, marker="o", markersize=5, label=result["label"],
        )  # fmt: skip
        prediction = result["many_qubits"]["prediction"]
        axes.plot(
            [prediction["n"]], [prediction["predicted_s"]],
            linestyle="none", marker="o", markersize=9, markerfacecolor="none",
            markeredgecolor=color, markeredgewidth=1.5,
        )  # fmt: skip
    handles, _ = axes.get_legend_handles_labels()
    handles.append(Line2D(
        [], [], linestyle="none", marker="o", markersize=9, markerfacecolor="none",
        markeredgecolor="#555555", label="predicted before the run",
    ))  # fmt: skip
    axes.legend(handles=handles, frameon=False)
    axes.set_yscale("log")
    axes.set_xlabel("qubits n")
    axes.set_ylabel("wall time (s)")
    axes.set_title(f"Brickwork circuit, depth {results[0]['depth']}, one noiseless shot")
    axes.xaxis.set_major_locator(MaxNLocator(integer=True))
    axes.grid(True, which="major", color="#e0e0e0", linewidth=0.8)
    axes.set_axisbelow(True)
    for side in ("top", "right"):
        axes.spines[side].set_visible(False)
    return figure


def summary(results: list[dict]) -> list[dict]:
    """Per version: the largest-n point, shots/s, and speed-ups against results[0]."""
    reference = results[0]
    reference_s = {p["n"]: p["wall_s"] for p in reference["many_qubits"]["points"]}
    rows = []
    for result in results:
        points = result["many_qubits"]["points"]
        largest = points[-1]
        common_n = max(p["n"] for p in points if p["n"] in reference_s)
        wall_s = next(p["wall_s"] for p in points if p["n"] == common_n)
        rows.append({
            "label": result["label"],
            "largest_n": largest["n"],
            "peak_rss_bytes": largest["peak_rss_bytes"],
            "effective_bandwidth_gbps": largest["effective_bandwidth_gbps"],
            "fraction_of_peak": largest["fraction_of_peak"],
            "shots_per_s": result["many_shots"]["shots_per_s"],
            "common_n": common_n,
            "many_qubits_speed_up": reference_s[common_n] / wall_s,
            "many_shots_speed_up": (
                result["many_shots"]["shots_per_s"] / reference["many_shots"]["shots_per_s"]
            ),
        })  # fmt: skip
    return rows


def plot(results_dir: Path) -> None:
    """Writes time_vs_n.png and prints the summary table."""
    results = load_results(results_dir)
    path = results_dir / "time_vs_n.png"
    figure = draw(results)
    figure.savefig(path, dpi=150)
    plt.close(figure)
    print(path)
    print(
        f"{'version':<16} {'max n':>5} {'peak RSS':>10} {'eff. BW':>12} {'% peak':>7} "
        f"{'shots/s':>10} {'speed-up (n)':>14} {'shots speed-up':>14}"
    )
    for row in summary(results):
        print(
            f"{row['label']:<16} {row['largest_n']:>5} "
            f"{row['peak_rss_bytes'] / 2**30:>6.3g} GiB "
            f"{row['effective_bandwidth_gbps']:>7.3g} GB/s "
            f"{100 * row['fraction_of_peak']:>6.1f}% "
            f"{row['shots_per_s']:>10.4g} "
            f"{row['many_qubits_speed_up']:>8.2f}x ({row['common_n']:>2}) "
            f"{row['many_shots_speed_up']:>13.2f}x"
        )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)

    run_parser = commands.add_parser("run", help="measure one version in both regimes")
    run_parser.add_argument("label", help="version name, e.g. reference")
    run_parser.add_argument("--sim", type=Path, required=True, help="sim executable")
    run_parser.add_argument("--commit", required=True, help="commit the sim was built from")
    run_parser.add_argument(
        "--max-n", type=int, required=True, help="largest qubit count (at least 5)"
    )
    run_parser.add_argument("--peak-bandwidth-gbps", type=float, required=True)
    run_parser.add_argument("--compiler", required=True, help="e.g. 'Apple clang 21'")
    run_parser.add_argument("--compiler-flags", required=True, help="e.g. '-O3 -DNDEBUG'")
    run_parser.add_argument("--results-dir", type=Path, default=RESULTS_DIR)

    plot_parser = commands.add_parser("plot", help="plot every result in the directory")
    plot_parser.add_argument("--results-dir", type=Path, default=RESULTS_DIR)

    args = parser.parse_args(argv)
    if args.command == "plot":
        plot(args.results_dir)
    elif args.command == "run":
        if args.max_n < 5:
            run_parser.error("--max-n must be at least 5 to predict the largest n")
        run(args)
    return 0


if __name__ == "__main__":
    sys.exit(main())
