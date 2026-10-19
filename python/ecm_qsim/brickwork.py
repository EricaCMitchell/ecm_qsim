"""Brickwork benchmark circuits: layers of random RY rotations and CNOT bricks.

Each layer is one time step of RY on every qubit, then one time step of CNOTs
on neighbouring pairs. Odd layers pair (0, 1), (2, 3), ...; even layers pair
(1, 2), (3, 4), ....
"""

import numpy as np


def brickwork_circuit(n: int, depth: int, seed: int) -> str:
    """Circuit text for n qubits and depth layers, with angles drawn from seed.

    The circuit has exactly 2 * depth time steps. A CNOT step with no pairs
    (n = 2, even layer) is kept as an empty step.
    """
    if n < 2:
        raise ValueError(f"need at least 2 qubits, got {n}")
    if depth < 1:
        raise ValueError(f"depth must be at least 1, got {depth}")

    rng = np.random.default_rng(seed)
    lines = [f"# brickwork n={n} depth={depth} seed={seed}", f"QUBITS {n}"]
    for layer in range(1, depth + 1):
        angles = rng.uniform(0, 2 * np.pi, size=n)
        # float() so that repr gives the plain number, not np.float64(...).
        lines += [f"RY {float(theta)!r} {q}" for q, theta in enumerate(angles)]
        lines.append("TICK")
        start = 0 if layer % 2 == 1 else 1
        lines += [f"CNOT {q} {q + 1}" for q in range(start, n - 1, 2)]
        lines.append("TICK")
    return "\n".join(lines) + "\n"
