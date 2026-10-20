# Fusion profile (fusion 35b3578, Apple M3, Apple clang 21, -O3 -DNDEBUG)

Time Profiler via `xcrun xctrace record`, 1 ms samples, inclusive share of samples. Same
procedure, circuits and options as `specialized-profile.md`.
Traces (local only): fusion-many-qubits.trace, fusion-many-shots.trace.

| Function group                             | n=29 noiseless, 1 shot (101918 samples) | n=12 noisy, 100000 shots (26707 samples) |
|--------------------------------------------|------:|------:|
| StateVector::apply_real_2q                 | 95.4% | 83.0% |
| StateVector::apply_real_1q                 | 1.9%  | 5.0%  |
| StateVector::apply_2q (complex clusters)   | 0%    | 0.2%  |
| fuse_trajectory (fused-matrix construction) | 0%   | 8.0%  |
| allocation (StateVector ctor, __bzero)     | 0.4%  | 0.6%  |
| measurement (probabilities + sample_index) | 2.3%  | 1.7%  |
| sample_pauli_errors                        | 0%    | 0.7%  |

Gate application takes about 97% of samples at n = 29 and 88% at n = 12. Fused-matrix
construction is negligible at n = 29 and takes 8.0% of the many-shots samples. Clang reports
that the `apply_real_2q` inner loop is not vectorized. The `apply_ry` and `apply_real_1q` inner
loops are vectorized (width 2).

## Results (fusion.json, run 2026-10-20)

Same machine, circuits, seeds and peak bandwidth (100 GB/s) as `reference.json` and
`specialized.json`. The reference is `v1-ref` (3cd8956). Correctness against it at n = 12 and
10000 shots: noiseless counts are an exact match and the noisy chi-square p-value is 1.0.

### Many qubits (noiseless, 1 shot)

Greedy fusion to width two turns depth-20 brickwork into `10 n - 1` fused gates, and each gate
is one full pass over the state: 32 logical bytes per amplitude. The fusion logical-byte model
is `(10 n - 1) * 32 * 2^n`. The `modeled_bytes` and bandwidth fields in `fusion.json` retain the
runner's per-gate convention, so the table derives the fusion values from the recorded wall
times.

| n | Wall time (s) | Speed-up vs specialized | Speed-up vs v1-ref | Logical-byte rate (GB/s) | % of peak |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 24 | 2.628 | 1.77x | 3.54x | 48.8 | 48.8 |
| 26 | 11.31 | 1.76x | 3.57x | 49.2 | 49.2 |
| 28 | 48.61 | 1.78x | 3.53x | 49.3 | 49.3 |
| 29 | 103.7 | 1.66x | 3.46x | 47.9 | 47.9 |

Peak memory at n = 29 is 12.01 GiB, unchanged. The runner predicted 100.7 s for n = 29 from the
n = 28 time; the measurement is 103.7 s (+2.9%). The traffic-model prediction written before
the implementation was 69.1 s (2.49x); the measurement is 50% slower (1.66x).

Fused real dense arithmetic is `10 (n - 1) * 14 + 9 * 6` FLOPs per amplitude, or 2.134e12 FLOPs
at n = 29: 20.6 GFLOP/s, against 10.8 GFLOP/s for `specialized`.

### Many shots (n = 12, 100000 shots, p_gate = p_idle = 1e-3, p_readout = 1e-2)

| Wall time | Shots/s | Time per shot | Peak memory | Speed-up vs specialized | Speed-up vs v1-ref |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 27.19 s | 3678 | 0.272 ms | 1.6 MiB | 1.07x | 3.77x |

Useful arithmetic is `10 (n - 1) * 14 + 9 * 6` = 1594 FLOPs per amplitude: 24.0 GFLOP/s, against
20.2 GFLOP/s for `specialized` (1440 FLOPs per amplitude).

At large n the fused version moves about 48-49 GB/s of logical bytes, less than the 69-72 GB/s
of the specialized kernels, while its FLOP rate doubled to the same 20-24 GFLOP/s range as the
cache-resident many-shots regime. Thus the non-vectorized 4x4 kernel, not DRAM bandwidth, now
limits both regimes.
