# ecm_qsim

A noisy quantum circuit simulator. It uses a state vector and Monte Carlo trajectories for noise.

## Requirements

- C++20 compiler (tested with Apple clang 21)
- CMake 3.24 or later, and Ninja
- Eigen 3.4 or later (`brew install eigen`)
- Catch2 3 (`brew install catch2`), optional: if it is not installed, CMake downloads it
- Python 3.11 or later, for tooling and tests only

## Build

```sh
cmake -S . -B build -G Ninja
cmake --build build
ctest --test-dir build
```

If Catch2 is not installed, the first configure downloads it and needs network access. The default
build type is `Release`.

## Python tooling

```sh
python3 -m venv .venv
source .venv/bin/activate
pip install -e './python[dev]'
pytest tests/python
```

## Usage

```sh
build/sim examples/ghz3.txt --shots 1000 --p-gate 0.001 --p-idle 0.001 --p-readout 0.01 --seed 1
```

`--shots` and `--seed` are required. The three noise probabilities default to 0. The output has one
line per measured bitstring, with qubit 0 as the rightmost character:

    000 496
    001 9
    ...

`--json` prints `{"counts": {"000": 496, ...}}` instead. The exit status is 0 on success, 1 if the
circuit cannot be read or run, and 2 for invalid arguments.

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
tests/python/       Python tests (pytest)
python/ecm_qsim/    Python tooling
benchmarks/         benchmark scripts and results
experiments/        noise validation experiments
examples/           example circuits
```
