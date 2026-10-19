"""Idle-noise experiment: fit p_idle from the decay of a qubit prepared in 1.

Runs the sim binary for each m, prints the data and the fit, and saves a
figure. p_gate is 0 and p_readout is treated as known.
"""

import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from ecm_qsim import idle, sim  # noqa: E402

DATA_COLOR = "#2a78d6"
FIT_COLOR = "#eb6834"
INJECTED_COLOR = "#52514e"
TEXT_COLOR = "#0b0b0b"
MUTED_COLOR = "#b5b4ae"


def plot(
    data: idle.IdleData,
    fit: idle.FitResult,
    p_idle: float,
    p_readout: float,
    out: Path,
) -> None:
    """Saves the data with the fitted and injected curves."""
    m_smooth = np.linspace(0, data.m.max(), 400)

    fig, ax = plt.subplots(figsize=(7, 4.5), layout="constrained")
    ax.axhline(0.5, color=MUTED_COLOR, linestyle=":", linewidth=1, label="limit 1/2")
    ax.plot(
        m_smooth,
        idle.analytic_p(m_smooth, p_idle, p_readout),
        color=INJECTED_COLOR,
        linestyle="--",
        linewidth=1.2,
        label=f"analytic, injected p_idle = {p_idle:g}",
    )
    ax.plot(
        m_smooth,
        idle.analytic_p(m_smooth, fit.p_idle, p_readout),
        color=FIT_COLOR,
        linewidth=1.5,
        label=f"fit, p_idle = {fit.p_idle:.4g} ± {fit.stderr:.2g}",
    )
    ax.errorbar(
        data.m,
        data.p,
        yerr=data.p_stderr,
        fmt="o",
        markersize=5,
        color=DATA_COLOR,
        elinewidth=1,
        label=f"simulated, {data.shots[0]} shots (±1 s.e.)",
    )
    ax.set_xlabel("idle steps m")
    ax.set_ylabel("fraction of shots returning 1")
    ax.set_title(
        f"Idle noise: p_readout = {p_readout:g} (known), "
        f"χ²/dof = {fit.chi2:.1f}/{fit.dof}",
        color=TEXT_COLOR,
    )
    ax.legend(frameon=False, labelcolor=TEXT_COLOR)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(color=MUTED_COLOR, linewidth=0.5, alpha=0.5)
    fig.savefig(out, dpi=150)
    plt.close(fig)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--p-idle", type=float, default=0.01)
    parser.add_argument("--p-readout", type=float, default=0.01)
    parser.add_argument("--shots", type=int, default=10_000, help="shots per m value")
    parser.add_argument("--seed", type=int, default=1, help="m value i uses seed + i")
    parser.add_argument(
        "--m", type=int, nargs="+", help="idle steps (default: 0 to 3 decay lengths)"
    )
    parser.add_argument("--sim", type=Path, default=sim.default_sim_path())
    parser.add_argument("--out", type=Path, default=Path("experiments/idle_noise.png"))
    args = parser.parse_args(argv)

    if args.m is None:
        try:
            m = idle.default_m_values(args.p_idle)
        except ValueError as error:
            parser.error(f"{error}; give --m")
    else:
        m = np.asarray(args.m)

    predicted = idle.fit_p_idle(
        m, idle.analytic_p(m, args.p_idle, args.p_readout), args.shots, args.p_readout
    )
    print(f"predicted stderr of p_idle for {args.shots} shots: {predicted.stderr:.3g}")

    try:
        data = idle.run(
            m,
            p_idle=args.p_idle,
            p_readout=args.p_readout,
            shots=args.shots,
            seed=args.seed,
            sim_path=args.sim,
        )
    except sim.SimError as error:
        print(error, file=sys.stderr)
        return 1
    fit = idle.fit_p_idle(data.m, data.p, data.shots, args.p_readout)

    print(f"{'m':>6} {'ones':>8} {'shots':>8}  p ± stderr")
    for row in zip(data.m, data.ones, data.shots, data.p, data.p_stderr):
        print(f"{row[0]:>6} {row[1]:>8} {row[2]:>8}  {row[3]:.4f} ± {row[4]:.4f}")
    print(f"p_idle fit = {fit.p_idle:.5g} ± {fit.stderr:.2g}")
    print(f"p_idle injected = {args.p_idle:g}")
    print(f"z = {(fit.p_idle - args.p_idle) / fit.stderr:+.2f}")
    print(f"chi2/dof = {fit.chi2:.1f}/{fit.dof}")

    plot(data, fit, args.p_idle, args.p_readout, args.out)
    print(f"figure: {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
