# Reference profile (v1-ref 3cd8956, Apple M3, Apple clang 21, -O3 -DNDEBUG)

Time Profiler via `xcrun xctrace record`, 1 ms samples, inclusive share of samples.
Traces: reference-many-qubits.trace, reference-many-shots.trace.

| Function group                         | n=29 noiseless, 1 shot (364182 samples) | n=12 noisy, 100000 shots (102692 samples) |
|----------------------------------------|------:|------:|
| StateVector::apply_1q                  | 59.8% | 63.5% |
| StateVector::apply_2q                  | 39.2% | 35.1% |
| matrix construction (matrix_1q/2q)     | <0.05% | 0.3% |
| allocation (StateVector ctor, malloc)  | 0.1% (__bzero) | 0.2% |
| measurement (probabilities + sample_index) | 0.8% | 0.5% |

Both regimes spend about 99% of samples in the two gate kernels. In the many-shots regime,
the state is rebuilt and the circuit re-run per shot, so the per-shot kernel work dominates.

## Reference results (reference.json, run 2026-10-20T02:39:50Z)

Same binary and machine as the profiles: Apple M3, 8 physical cores, 24 GiB, macOS 26.6.2,
no GPU. Depth 20, angle seed 42, sim seed 42, one run for each point. Peak memory is the peak
physical footprint (compressed pages included). Peak bandwidth supplied: 100 GB/s.

### Many qubits (noiseless, 1 shot)

Modeled bytes: one read and one write of every 16-byte amplitude for each gate,
`(20 n + 10 (n - 1)) * 32 * 2^n`. Effective bandwidth = modeled bytes / wall time.

| n | Wall time (s) | Peak memory | Modeled bytes | Effective bandwidth (GB/s) | % of peak |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 4 | 0.004052 | 0.9 MiB | 5.632e+04 | 0.0139 | 0.0 |
| 5 | 0.003015 | 1.0 MiB | 1.434e+05 | 0.0475 | 0.0 |
| 6 | 0.003059 | 1.0 MiB | 3.482e+05 | 0.114 | 0.1 |
| 7 | 0.002883 | 1.0 MiB | 8.192e+05 | 0.284 | 0.3 |
| 8 | 0.002431 | 1.0 MiB | 1.884e+06 | 0.775 | 0.8 |
| 9 | 0.002588 | 1.0 MiB | 4.260e+06 | 1.65 | 1.6 |
| 10 | 0.002723 | 1.0 MiB | 9.503e+06 | 3.49 | 3.5 |
| 11 | 0.002902 | 1.0 MiB | 2.097e+07 | 7.23 | 7.2 |
| 12 | 0.003314 | 1.1 MiB | 4.588e+07 | 13.8 | 13.8 |
| 13 | 0.004378 | 1.2 MiB | 9.961e+07 | 22.8 | 22.8 |
| 14 | 0.007332 | 1.4 MiB | 2.150e+08 | 29.3 | 29.3 |
| 15 | 0.01369 | 1.8 MiB | 4.614e+08 | 33.7 | 33.7 |
| 16 | 0.02668 | 2.5 MiB | 9.857e+08 | 36.9 | 36.9 |
| 17 | 0.05152 | 4.0 MiB | 2.097e+09 | 40.7 | 40.7 |
| 18 | 0.1068 | 7.1 MiB | 4.446e+09 | 41.6 | 41.6 |
| 19 | 0.2294 | 13.0 MiB | 9.395e+09 | 41.0 | 41.0 |
| 20 | 0.4779 | 25.1 MiB | 1.980e+10 | 41.4 | 41.4 |
| 21 | 1.010 | 49.1 MiB | 4.161e+10 | 41.2 | 41.2 |
| 22 | 2.154 | 97.1 MiB | 8.724e+10 | 40.5 | 40.5 |
| 23 | 4.475 | 193.1 MiB | 1.825e+11 | 40.8 | 40.8 |
| 24 | 9.287 | 385.2 MiB | 3.812e+11 | 41.0 | 41.0 |
| 25 | 19.64 | 769.4 MiB | 7.946e+11 | 40.5 | 40.5 |
| 26 | 40.35 | 1.50 GiB | 1.654e+12 | 41.0 | 41.0 |
| 27 | 83.82 | 3.00 GiB | 3.436e+12 | 41.0 | 41.0 |
| 28 | 171.4 | 6.00 GiB | 7.130e+12 | 41.6 | 41.6 |
| 29 | 359.0 | 12.01 GiB | 1.477e+13 | 41.2 | 41.2 |

Prediction for n = 29, printed before the run from the n = 28 time: 355.3 s. Measured: 359.0 s
(+1.1 %).

The algorithm requires 1.868e12 FLOPs for RY (6 per amplitude) and none for CNOT. The generic
reference kernels actually execute dense complex matrix-vector products: about 14 real FLOPs per
amplitude for RY and 30 for CNOT, or 8.869e12 FLOPs for this circuit. This distinction matters for
the optimized-kernel plans, which remove work represented by zero matrix entries.

### Many shots (n = 12, 100000 shots, p_gate = p_idle = 1e-3, p_readout = 1e-2)

| Wall time | Shots/s | Time per shot | Peak memory |
| ---: | ---: | ---: | ---: |
| 102.5 s | 975.8 | 1.025 ms | 1.5 MiB |

Derived: with the n = 12 byte model (4.588e7 bytes per noiseless circuit), 100000 shots move
4.6e12 modeled bytes, an effective 44.8 GB/s. The 64 KiB state fits in cache, so this figure is
not a DRAM bandwidth. The model excludes the injected Pauli errors: 480 qubit-steps per shot at
p = 1e-3 give about 0.5 extra gates per shot, against 350 circuit gates.
