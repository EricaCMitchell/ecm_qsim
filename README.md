# ecm_qsim

A noisy quantum circuit simulator. It uses a state vector and Monte Carlo trajectories for noise.

## Requirements

- C++20 compiler (tested with Apple clang 21)
- CMake 3.24 or later, and Ninja
- Eigen 3.4 or later (`brew install eigen`)
- Catch2 3 (`brew install catch2`), optional: if it is not installed, CMake downloads it

## Build

```sh
cmake -S . -B build -G Ninja
cmake --build build
ctest --test-dir build
```

If Catch2 is not installed, the first configure downloads it and needs network access. The default
build type is `Release`.

## Usage

```sh
build/sim examples/ghz3.txt --shots 1000 --p-gate 0.001 --p-idle 0.001 --p-readout 0.01 --seed 1
```

`--shots` and `--seed` are required. The three noise probabilities default to 0. The output has one
line per observed bitstring, with qubit 0 as the rightmost character:

    000 496
    001 9
    ...

`--json` prints `{"counts": {"000": 496, ...}}` instead. The exit status is 0 on success, 1 if the
circuit cannot be read or run, and 2 for invalid arguments.

`--p-gate` and `--p-idle` are the probabilities of a uniformly selected Pauli X, Y or Z error on
each active or idle qubit after a time step. `--p-readout` is the independent probability of
flipping each measured bit. All probabilities must be in `[0, 1]`.

## Circuit format

The first non-comment statement declares the number of qubits:

```text
QUBITS 3
```

The supported one-qubit gates are `X`, `Y`, `Z`, `H`, `S` and `T`:

```text
H 0
X 1
```

Rotation angles are in radians:

```text
RX 1.5707963267948966 0
RY -0.5 1
RZ 0.785 2
```

The supported two-qubit gates are `CNOT`, `CZ` and `CY`. The first qubit is the control and the
second is the target:

```text
CNOT 0 1
```

Operations remain in the same time step until `TICK`; a qubit may occur in at most one operation
per time step. Blank lines and lines whose first token begins with `#` are ignored. Inline comments
are not supported. The final time step is closed at end of file, so a trailing `TICK` is optional.
See [`examples/ghz3.txt`](examples/ghz3.txt) for a complete circuit.

## Reproducibility

A seed gives the same random draws on every platform. The generator is `std::mt19937_64`, and the
simulator uses no `std::` distributions, whose output depends on the standard library. Counts can
differ between compilers only when a draw falls within rounding error of the boundary between two
outcomes.

## Layout

```text
include/ecm_qsim/   public headers
src/                library sources and the sim CLI (main.cpp)
tests/cpp/          C++ unit tests (Catch2)
tests/python/       reserved for future Python tests (currently empty)
python/ecm_qsim/    Python package scaffold (no tooling implemented yet)
benchmarks/         reserved benchmark-results directory
experiments/        reserved experiments directory
examples/           example circuits
```

Python packaging metadata exists under `python/`, but the repository does not currently contain
Python utilities or Python tests.
