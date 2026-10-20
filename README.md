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

## Python tooling

The Python package in `python/` runs `build/sim` as a subprocess. Install it in a virtual
environment at the repository root:

```sh
python3 -m venv .venv
.venv/bin/pip install -e './python[dev]'
.venv/bin/pytest tests/python
```

Set `ECM_QSIM_SIM` to use a `sim` binary at a different path. Tests that need the binary are
skipped if it is not built.

### Idle-noise experiment

```sh
.venv/bin/python experiments/idle_noise.py --p-idle 0.01 --p-readout 0.01 --shots 10000
```

The script prepares one qubit in 1, leaves it idle for `m` time steps with `p_gate = 0`, and
measures it. It fits `p_idle` to the fraction of shots that return 1, with `p_readout` known:

    p(m) = 1/2 + (1 - 2 p_readout)/2 * (1 - 4 p_idle / 3)^m

An X or Y error flips the bit and a Z error does not, so each idle step flips the bit with
probability `2 p_idle / 3`. The script prints the data, the fitted `p_idle` with its standard
error, and the stderr predicted before the run. It saves the figure to
`experiments/idle_noise.png` (`--out` to change). Use `--m` to give the idle step counts.

The figures in `experiments/figures/` come from these commands:

```sh
P=.venv/bin/python
$P experiments/idle_noise.py --p-idle 0.01 --out experiments/figures/idle_noise_p0.01.png
$P experiments/idle_noise.py --p-idle 0.9 --m $(seq -s ' ' 0 10) \
    --out experiments/figures/idle_noise_p0.9.png
```

### Brickwork benchmark circuits

```sh
.venv/bin/python benchmarks/brickwork_generator.py --qubits 12 --depth 20
```

Each of the `--depth` layers is one time step of `RY` with a random angle on every qubit, then
one time step of CNOTs on neighbouring pairs: `(0, 1), (2, 3), …` in odd layers and
`(1, 2), (3, 4), …` in even layers. The script writes
`benchmarks/circuits/brickwork_n<n>_d<D>.txt` (`--out` to change) and prints the path. These
files are not committed. Each run draws new angles unless `--seed` is given; the seed is written
in the file header.

### Benchmarks

`benchmarks/run_benchmarks.py` measures a `sim` binary in two regimes and plots the results. The
circuits are brickwork circuits of depth 20 with angle seed 42, and the simulator seed is 42:

1. Many qubits: noiseless, one shot, `n = 4` to `--max-n`.
2. Many shots: `n = 12`, `p_gate = p_idle = 0.001`, `p_readout = 0.01`, 100000 shots.

The script runs on macOS and Linux. It needs CMake, a C++ compiler, and Python 3.11 or later
with the `.venv` from [Python tooling](#python-tooling). On Linux, it also needs `lscpu`.

Build the reference version (tag `v1-ref`) as Release in a separate directory:

```sh
mkdir -p build/reference/src
git archive v1-ref | tar -x -C build/reference/src
cmake -S build/reference/src -B build/reference/build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build/reference/build --target sim
```

Build each optimized version as Release in its own directory, for example
`cmake -S . -B build/release -G Ninja -DCMAKE_BUILD_TYPE=Release`. Use the same compiler and
optimization flags for all versions. CMake records them in `CMakeCache.txt`
(`CMAKE_CXX_COMPILER` and `CMAKE_CXX_FLAGS_RELEASE`).

Measure the reference, then plot:

```sh
P=.venv/bin/python
$P benchmarks/run_benchmarks.py run reference --sim build/reference/build/sim \
    --commit "$(git rev-parse 'v1-ref^{commit}')" --max-n 29 --peak-bandwidth-gbps 100 \
    --compiler "Apple clang 21.0.0 (clang-2100.3.34.2)" --compiler-flags "-O3 -DNDEBUG"
$P benchmarks/run_benchmarks.py plot
```

`--max-n` is the largest `n` that fits in memory. The state vector needs `16 * 2^n` bytes, so
`n = 29` needs 8 GiB. `--peak-bandwidth-gbps` is the theoretical peak memory bandwidth of the
machine (100 GB/s for an Apple M3). Before the run at `--max-n`, the script prints a prediction
of its time from the time at `--max-n - 1`. If a run fails, the script stops and writes no result.

To measure an optimized version, give it its own label and binary:

```sh
$P benchmarks/run_benchmarks.py run simd --sim build/release/sim \
    --reference-sim build/reference/build/sim \
    --commit "$(git rev-parse HEAD)" --max-n 29 --peak-bandwidth-gbps 100 \
    --compiler "Apple clang 21.0.0 (clang-2100.3.34.2)" --compiler-flags "-O3 -DNDEBUG"
$P benchmarks/run_benchmarks.py plot
```

Before timing an optimized binary, `run` checks it against the reference binary at `n = 12`
with 10000 shots. The noiseless counts must match exactly. For the noisy counts, `run` puts each
outcome in a group by its noiseless count (0, 1, 2, 3, or 4 and more). A two-sample chi-square
test over these groups must give `p >= 0.001`. Noise moves shots from likely to unlikely
outcomes, and the groups make this move visible. `run` also rejects counts that do not add up to
the requested shots.

`run` writes `benchmarks/results/<LABEL>.json`. `plot` reads all results in that folder and
writes `benchmarks/results/time_vs_n.png`. It also prints a table with peak memory, effective
bandwidth, shots per second, and the speed-up against the reference. Curated JSON, PNG, and
Markdown evidence is committed; raw Instruments and `perf` profiles are not.

On macOS, peak memory is the peak physical footprint of `sim`. This value includes pages that
macOS compresses. The maximum resident set size does not include them, so it is too low near the
memory limit. On Linux, peak memory is the maximum resident set size. If Linux swaps pages out,
this value is less than the memory that `sim` allocated.

#### Profiling

Profile on the machine that made the results, with the same Release binary. Make the circuit
with `benchmarks/brickwork_generator.py --qubits 29 --depth 20 --seed 42`. On macOS, install
Xcode and select it with `sudo xcode-select -s /Applications/Xcode.app`. Then record a Time
Profiler trace and open it in Instruments:

```sh
xcrun xctrace record --template 'Time Profiler' \
    --output benchmarks/results/reference-many-qubits.trace \
    --launch -- build/reference/build/sim benchmarks/circuits/brickwork_n29_d20.txt \
    --shots 1 --seed 42
```

On Linux, you must have permission to run `perf` (see `/proc/sys/kernel/perf_event_paranoid`).
Record and read the profile:

```sh
perf record -o benchmarks/results/reference-many-qubits.perf.data -- \
    build/reference/build/sim benchmarks/circuits/brickwork_n29_d20.txt --shots 1 --seed 42
perf report --stdio -i benchmarks/results/reference-many-qubits.perf.data
```

For the many-shots profile, use the `n = 12` circuit, the name `reference-many-shots`, and the
options `--shots 100000 --seed 42 --p-gate 0.001 --p-idle 0.001 --p-readout 0.01`.

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
tests/python/       Python tests (pytest)
python/ecm_qsim/    Python tooling: sim subprocess wrapper, idle-noise model and fit,
                    brickwork circuits
benchmarks/         brickwork circuit generator and benchmark runner
experiments/        experiment scripts and figures
examples/           example circuits
```
