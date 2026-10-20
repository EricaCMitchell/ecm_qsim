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

The two specialized kernels take about 98% of samples in the n = 29 regime and 96% in the
n = 12 regime. `apply_ry` is now the larger share in both. The general kernels remain only for
the sampled Pauli errors.

## Results (specialized.json, run 2026-10-20)

Same machine, circuits, seeds and peak bandwidth (100 GB/s) as `reference.json`: Apple M3,
8 physical cores, 24 GiB, macOS 26.6.2, no GPU. The reference is `v1-ref` (3cd8956).
Correctness against it at n = 12 and 10000 shots: noiseless counts are an exact match and the
noisy chi-square p-value is 1.0.

### Many qubits (noiseless, 1 shot)

RY reads and writes every 16-byte amplitude, or 32 logical bytes per amplitude. The specialized
CNOT swaps only half of the amplitudes, or 16 logical bytes per amplitude in the full state. At
depth 20 the circuit has `20 n` RY gates and `10 (n - 1)` CNOT gates, so the specialized
logical-byte model is `(20 n * 32 + 10 (n - 1) * 16) * 2^n`. Effective bandwidth is modeled
bytes divided by wall time.

| n | Wall time (s) | Speed-up | Peak memory | Modeled bytes | Effective bandwidth (GB/s) | % of peak |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 4 | 0.003770 | 1.07x | 0.9 MiB | 4.864e+04 | 0.0129 | 0.0 |
| 5 | 0.003525 | 0.86x | 1.0 MiB | 1.229e+05 | 0.0349 | 0.0 |
| 6 | 0.002979 | 1.03x | 1.0 MiB | 2.970e+05 | 0.0997 | 0.1 |
| 7 | 0.003212 | 0.90x | 1.0 MiB | 6.963e+05 | 0.217 | 0.2 |
| 8 | 0.003137 | 0.78x | 1.0 MiB | 1.597e+06 | 0.509 | 0.5 |
| 9 | 0.003116 | 0.83x | 1.0 MiB | 3.604e+06 | 1.16 | 1.2 |
| 10 | 0.003006 | 0.91x | 1.0 MiB | 8.028e+06 | 2.67 | 2.7 |
| 11 | 0.003109 | 0.93x | 1.0 MiB | 1.769e+07 | 5.69 | 5.7 |
| 12 | 0.003365 | 0.98x | 1.1 MiB | 3.867e+07 | 11.5 | 11.5 |
| 13 | 0.003577 | 1.22x | 1.2 MiB | 8.389e+07 | 23.4 | 23.4 |
| 14 | 0.004706 | 1.56x | 1.4 MiB | 1.809e+08 | 38.4 | 38.4 |
| 15 | 0.006622 | 2.07x | 1.8 MiB | 3.880e+08 | 58.6 | 58.6 |
| 16 | 0.01068 | 2.50x | 2.5 MiB | 8.284e+08 | 77.5 | 77.5 |
| 17 | 0.02023 | 2.55x | 4.0 MiB | 1.762e+09 | 87.1 | 87.1 |
| 18 | 0.03801 | 2.81x | 7.1 MiB | 3.733e+09 | 98.2 | 98.2 |
| 19 | 0.07534 | 3.04x | 13.0 MiB | 7.885e+09 | 105 | 104.7 |
| 20 | 0.1602 | 2.98x | 25.1 MiB | 1.661e+10 | 104 | 103.7 |
| 21 | 0.4006 | 2.52x | 49.1 MiB | 3.490e+10 | 87.1 | 87.1 |
| 22 | 0.9530 | 2.26x | 97.1 MiB | 7.315e+10 | 76.8 | 76.8 |
| 23 | 2.237 | 2.00x | 193.2 MiB | 1.530e+11 | 68.4 | 68.4 |
| 24 | 4.650 | 2.00x | 385.3 MiB | 3.194e+11 | 68.7 | 68.7 |
| 25 | 9.831 | 2.00x | 769.4 MiB | 6.657e+11 | 67.7 | 67.7 |
| 26 | 19.95 | 2.02x | 1.50 GiB | 1.385e+12 | 69.4 | 69.4 |
| 27 | 43.74 | 1.92x | 3.00 GiB | 2.878e+12 | 65.8 | 65.8 |
| 28 | 86.43 | 1.98x | 6.00 GiB | 5.970e+12 | 69.1 | 69.1 |
| 29 | 172.2 | 2.08x | 12.01 GiB | 1.237e+13 | 71.8 | 71.8 |

The `modeled_bytes` and bandwidth fields in `specialized.json` retain the benchmark runner's
original full-pass convention. The table above derives the specialized logical traffic from the
recorded wall times. At n = 29, RY performs 1.868e12 real FLOPs and CNOT performs none.

Prediction for n = 29, printed before the run from the n = 28 time: 179.1 s. Measured: 172.2 s
(-4%). The measured value is the median of three standalone reruns: 172.24 s, 174.06 s and
167.25 s. The runner's own n = 29 run took 253.0 s and remains in `specialized.json` as
`dropped_wall_s`. The three reruns agree within 4%, and the n = 29 profile took about 175 s.

### Many shots (n = 12, 100000 shots, p_gate = p_idle = 1e-3, p_readout = 1e-2)

| Wall time | Shots/s | Time per shot | Peak memory | Speed-up |
| ---: | ---: | ---: | ---: | ---: |
| 29.18 s | 3427 | 0.292 ms | 1.5 MiB | 3.51x |

The noiseless logical-byte model gives 3.867e12 bytes over 100000 shots, or 133 GB/s. This is not
a DRAM bandwidth: the 64 KiB state fits in the cache hierarchy. It excludes the approximately
0.5 sampled Pauli gates per shot.

At large n, the specialized rate is about 66-72% of the theoretical machine-wide bandwidth from
one thread; no hardware counters were collected to separate memory traffic from loop overhead.
Rates above 100 GB/s at n = 19 and 20 reflect cache reuse rather than DRAM traffic. Gate
application still takes about 96-98% of profile samples. With more time, the next optimization would
be gate fusion, reducing the number of full-state passes.
