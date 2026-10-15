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
