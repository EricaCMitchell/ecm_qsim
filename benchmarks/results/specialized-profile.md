# Specialized-kernel profile (specialized 794ea38, Apple M3, Apple clang 21, -O3 -DNDEBUG)

Time Profiler via `xcrun xctrace record`, 1 ms samples, inclusive share of samples. Same
procedure, circuits and options as `reference-profile.md`.
Traces (local only): specialized-many-qubits.trace, specialized-many-shots.trace.

| Function group                         | n=29 noiseless, 1 shot (175109 samples) | n=12 noisy, 100000 shots (29189 samples) |
|----------------------------------------|------:|------:|
| StateVector::apply_ry                  | 71.7% | 79.0% |
| StateVector::apply_cnot                | 26.3% | 17.1% |
| StateVector::apply_1q (Pauli errors)   | 0%    | 0.4% |
| StateVector::apply_2q                  | 0%    | 0%   |
| allocation (StateVector ctor, __bzero) | 0.3%  | 0.5% |
| measurement (probabilities + sample_index) | 1.7% | 1.5% |
| sample_pauli_errors                    | 0%    | 0.7% |

The two specialized kernels take about 98 % of samples in the n = 29 regime and 96 % in the
n = 12 regime. apply_ry is now the larger share in both. The general kernels remain only for the
sampled Pauli errors.

## Results (specialized.json, run 2026-10-20)

Same machine, circuits, seeds and peak bandwidth (100 GB/s) as reference.json. The reference sim
is v1-ref (3cd8956). Correctness against the reference at n = 12, 10000 shots: noiseless counts
are an exact match, and the noisy chi-square p-value is 1.0. The kernels do not change the draw
order, so the noisy counts are also identical.

| n | Reference (s) | Specialized (s) | Speed-up | Peak memory |
| ---: | ---: | ---: | ---: | ---: |
| 12 | 0.003314 | 0.003365 | 0.98x | 1.1 MiB |
| 16 | 0.02668 | 0.01068 | 2.50x | 2.5 MiB |
| 20 | 0.4779 | 0.1602 | 2.98x | 25.1 MiB |
| 24 | 9.287 | 4.650 | 2.00x | 385.3 MiB |
| 26 | 40.35 | 19.95 | 2.02x | 1.50 GiB |
| 27 | 83.82 | 43.74 | 1.92x | 3.00 GiB |
| 28 | 171.4 | 86.43 | 1.98x | 6.00 GiB |
| 29 | 359.0 | 172.2 | 2.08x | 12.01 GiB |

Prediction for n = 29, printed before the run from the n = 28 time: 179.1 s. Measured: 172.2 s
(-4 %). The n = 29 value is the median of three standalone reruns (`/usr/bin/time -l`, same binary
and arguments): 172.24 s, 174.06 s and 167.25 s. The runner's own n = 29 run took 253.0 s (+41 %),
and specialized.json keeps it as `dropped_wall_s`. The three reruns agree within 4 %, and the n = 29
profile (175109 samples, about 175 s) agrees with them, so that run was an outlier. Its cause is
not known: the reruns had no swaps and the same 12.9 GB peak footprint.

Many shots (n = 12, 100000 shots): 29.18 s, 3427 shots/s (0.292 ms per shot), 3.51x the
reference. Peak memory 1.5 MiB.
