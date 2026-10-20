# Memory results (memory c6678cd, Apple M3, Apple clang 21, -O3 -DNDEBUG)

`run()` samples the state in place: `sample_index(const StateVector&, double)` computes each
`|amplitude|^2` as it reads it. Before, `run()` sampled a `probabilities()` copy of 8 bytes per
amplitude. This change targets capacity, not speed. The kernels and fusion do not change, and
measurement is 2.3% of the n = 29 profile (`fusion-profile.md`). Thus there is no new profile.

## Prediction

Peak memory falls from `24 * 2^n` bytes (amplitudes and copy) to `16 * 2^n` bytes (amplitudes
only): from about 12 GiB to 8 GiB at n = 29. Wall time stays within run-to-run noise.

## Results (memory.json, run 2026-10-20)

Same machine, circuits, seeds and peak bandwidth (100 GB/s) as `fusion.json`. The reference is
`v1-ref` (3cd8956).

### Many qubits (noiseless, 1 shot)

| n | fusion peak (GiB) | `24 * 2^n` (GiB) | memory peak (GiB) | `16 * 2^n` (GiB) | fusion time (s) | memory time (s) |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 24 | 0.376 | 0.375 | 0.251 | 0.250 | 2.628 | 2.603 |
| 25 | 0.751 | 0.750 | 0.501 | 0.500 | 5.449 | 5.426 |
| 26 | 1.502 | 1.500 | 1.002 | 1.000 | 11.31 | 11.27 |
| 27 | 3.003 | 3.000 | 2.002 | 2.000 | 23.41 | 23.42 |
| 28 | 6.004 | 6.000 | 4.003 | 4.000 | 48.61 | 48.46 |
| 29 | 12.007 | 12.000 | 8.005 | 8.000 | 103.7 | 101.5 |

Both versions follow their models within 0.6% at n = 24 and within 0.2% from n = 26. At n = 29
the peak falls from 12.01 GiB to 8.01 GiB (-33.3%). The wall time is 101.5 s against 103.7 s
(-2.0%). This is a regression check, not a speed-up: the change is within run-to-run noise.

The qubit limit on this machine stays at 29. `StateVector::max_qubits()` now budgets 16 bytes
per amplitude in half of physical memory (12 GiB): 8 GiB fits at n = 29, and 16 GiB at n = 30
does not.

### Many shots (n = 12, 100000 shots, p_gate = p_idle = 1e-3, p_readout = 1e-2)

| Version | Wall time | Shots/s | Peak memory |
| --- | ---: | ---: | ---: |
| fusion | 27.19 s | 3678 | 1.63 MiB |
| memory | 26.65 s | 3752 | 1.67 MiB |

At n = 12 the copy is 32 KiB, so process overhead sets the peak. Shots/s changes by +2.0%,
which is within run-to-run noise.

## Correctness

- Against `v1-ref` at n = 12 and 10000 shots: the noiseless counts are an exact match, and the
  noisy chi-square p-value is 1.0.
- "run reference engine follows the documented draw order" still samples a `probabilities()`
  copy. It passes, so `run()` picks the same index as the old path for each of its 500 noisy
  shots.
- "sample_index on a state agrees with its probabilities" compares the two overloads just
  below, at and just above each cumulative boundary, for 1-qubit and 6-qubit random states.
