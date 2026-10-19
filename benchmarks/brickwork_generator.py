"""Write one brickwork benchmark circuit for a given number of qubits and layers.

Each layer is RY on every qubit with a random angle, then CNOTs on alternating
neighbouring pairs. Without --seed, each run draws new angles; the seed is
always written in the file header.
"""

import argparse
import secrets
import sys
from pathlib import Path

from ecm_qsim import brickwork

CIRCUITS_DIR = Path(__file__).resolve().parent / "circuits"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--qubits", type=int, required=True)
    parser.add_argument("--depth", type=int, default=20, help="number of layers")
    parser.add_argument("--seed", type=int, help="angle seed (default: random)")
    parser.add_argument(
        "--out",
        type=Path,
        help="output file (default: benchmarks/circuits/brickwork_n<n>_d<D>.txt)",
    )
    args = parser.parse_args(argv)

    seed = secrets.randbits(63) if args.seed is None else args.seed
    try:
        text = brickwork.brickwork_circuit(args.qubits, args.depth, seed)
    except ValueError as error:
        parser.error(str(error))

    out = args.out or CIRCUITS_DIR / f"brickwork_n{args.qubits}_d{args.depth}.txt"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(text)
    print(out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
