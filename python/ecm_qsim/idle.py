"""Idle-noise experiment: prepare 1, idle for m time steps, measure.

The fit recovers p_idle from the probability of reading 1 against m, with
p_readout known.
"""

from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from numpy.typing import ArrayLike
from scipy.optimize import curve_fit

from ecm_qsim import sim


def idle_circuit(m: int) -> str:
    """One qubit: X in the first time step, then m empty time steps."""
    if m < 0:
        raise ValueError(f"m must be non-negative, got {m}")
    return "QUBITS 1\nX 0\nTICK\n" + "TICK\n" * m


def analytic_p(m: ArrayLike, p_idle: float, p_readout: float) -> np.ndarray:
    """Analytic probability of reading 1 after m idle steps.

    p = 1/2 + (1 - 2 p_readout)/2 * (1 - 4 p_idle / 3)^m

    Per step, X or Y (2/3 of idle errors) flips the bit and Z does not, so
    the bit flips with probability 2 p_idle / 3.
    """
    decay = (1 - 4 * p_idle / 3) ** np.asarray(m, dtype=float)
    return 0.5 + (1 - 2 * p_readout) / 2 * decay


def binomial_stderr(p: ArrayLike, shots: ArrayLike) -> np.ndarray:
    """Binomial standard error of p, with p clipped to [1/n, 1 - 1/n].

    The clip keeps the error above 0 when every shot gives the same result.
    """
    shots = np.asarray(shots, dtype=float)
    p = np.clip(p, 1 / shots, 1 - 1 / shots)
    return np.sqrt(p * (1 - p) / shots)


def default_m_values(p_idle: float, points: int = 21) -> np.ndarray:
    """Evenly spaced integer m from 0 to three decay lengths."""
    if not 0 < p_idle < 0.75:
        raise ValueError(f"default m values need 0 < p_idle < 3/4, got {p_idle}")
    decay_length = -1 / np.log(1 - 4 * p_idle / 3)
    return np.unique(np.round(np.linspace(0, 3 * decay_length, points)).astype(int))


@dataclass(frozen=True)
class IdleData:
    m: np.ndarray
    shots: np.ndarray
    ones: np.ndarray

    @property
    def p(self) -> np.ndarray:
        """Simulated probability of reading 1."""
        return self.ones / self.shots

    @property
    def p_stderr(self) -> np.ndarray:
        return binomial_stderr(self.p, self.shots)


def run(
    m: Sequence[int],
    *,
    p_idle: float,
    p_readout: float,
    shots: int,
    seed: int,
    sim_path: Path | None = None,
) -> IdleData:
    """Runs one idle circuit per m with p_gate = 0. Point i uses seed + i."""
    ones = []
    for i, steps in enumerate(m):
        counts = sim.run(
            idle_circuit(int(steps)),
            shots=shots,
            seed=seed + i,
            p_idle=p_idle,
            p_readout=p_readout,
            sim_path=sim_path,
        )
        ones.append(counts.get("1", 0))
    return IdleData(
        m=np.asarray(m), shots=np.full(len(ones), shots), ones=np.asarray(ones)
    )


@dataclass(frozen=True)
class FitResult:
    p_idle: float
    stderr: float
    chi2: float
    dof: int


def fit_p_idle(
    m: ArrayLike, p: ArrayLike, shots: ArrayLike, p_readout: float
) -> FitResult:
    """Weighted least-squares fit of analytic_p to p, with p_readout known."""
    m = np.asarray(m)
    p = np.asarray(p, dtype=float)
    if not np.any(m > 0):
        raise ValueError("at least one m must be greater than 0")
    sigma = np.broadcast_to(binomial_stderr(p, shots), p.shape)
    # A large start makes the curve flat at every m > 0, so the fit cannot
    # move. Above 3/4 the curve alternates about 1/2, and only odd m tell
    # p_idle from 3/2 - p_idle.
    popt, pcov = curve_fit(
        lambda m, p_idle: analytic_p(m, p_idle, p_readout),
        m,
        p,
        p0=[1e-6],
        sigma=sigma,
        absolute_sigma=True,
        bounds=(0, 1),
    )
    p_idle = float(popt[0])
    chi2 = float(np.sum(((p - analytic_p(m, p_idle, p_readout)) / sigma) ** 2))
    return FitResult(p_idle, float(np.sqrt(pcov[0, 0])), chi2, len(m) - 1)
