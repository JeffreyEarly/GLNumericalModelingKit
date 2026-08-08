# Bundled versus native FFTW comparison

Benchmark status: **PASSED**
Horizontal bundled-FFTW readiness: **READY**

Native FFTW is a controlled reference only and is not a runtime requirement.

## Environment

| Field | Value |
|---|---|
| Run | 20260808T153906007Z-maca64-r2026a |
| MATLAB | 26.1.0.3312084 (R2026a) Update 4 |
| Architecture | maca64 |
| Operating system | Darwin 25.5.0 Darwin Kernel Version 25.5.0: Tue Jun  9 22:28:34 PDT 2026; root:xnu-12377.121.10~1/RELEASE_ARM64_T6050 arm64 |
| Processor | Apple M5 Max |
| Machine | MacBook Pro Mac17,6 |
| Physical memory | 48 GB |
| Threads | 18 |

## Engines

| Engine | Version | Resolved library | Matched alignment | Mismatch rejected | Unaligned |
|---|---|---|:---:|:---:|:---:|
| bundled-fftw | fftw-3.3.8 | /Applications/MATLAB_R2026a.app/bin/maca64/libmwfftw3.3.dylib | yes | yes | yes |
| native-fftw | fftw-3.3.11 | /Users/jearly/Documents/OceanKitRepositories/GLNumericalModelingKit/Matlab/Spectral/FFTW/build/native-fftw-3.3.11/install/lib/libfftw3.3.dylib | yes | yes | yes |

## Issue #38 discrepancy

Issue #38 selected one candidate from a large matrix using three screening samples, recreated that plan without an explicit wisdom reset, changed from a three-operation to a five-operation schedule, and did not persist the bundled finalist timings. The rows below replay the reconstructable finalists with isolated wisdom and full sampling; staged results are diagnostic only.

| Size | Role | Engine | Strategy | Planner | Alignment | #38 screen (ms) | Full replay (ms) | New matched result (ms) |
|---|---|---|---|---|---|---:|---:|---:|
| 256 x 256 x 64 | best-bundled-screening | bundled-fftw | guru-rank2 | measure | matched | 0.908 | 0.818 | 0.842 |
| 256 x 256 x 64 | selected-native-winner | native-fftw | guru-rank2 | exhaustive | unaligned | 0.827 | 0.729 | 0.692 |
| 1024 x 1024 x 30 | best-bundled-screening | bundled-fftw | staged-r2c-c2c | estimate | matched | 9.435 | 10.176 | 9.645 |
| 1024 x 1024 x 30 | selected-native-winner | native-fftw | guru-rank2 | patient | unaligned | 7.979 | 8.624 | 8.537 |

## Matched configuration comparison

`Bundled vs native total` is the percentage by which bundled time differs from native time; positive means bundled FFTW was slower.

| Size | Planner | Alignment | Status | MATLAB forward (ms) | Bundled raw (ms) | Bundled total (ms) | Native total (ms) | Bundled vs native total |
|---|---|---|---|---:|---:|---:|---:|---:|
| 256 x 256 x 64 | estimate | matched | passed | 1.895 | 1.087 | 1.144 | 1.056 | +8.3% |
| 256 x 256 x 64 | estimate | unaligned | passed | 2.033 | 1.070 | 1.143 | 1.091 | +4.7% |
| 256 x 256 x 64 | measure | matched | passed | 1.839 | 0.789 | 0.856 | 0.692 | +23.7% |
| 256 x 256 x 64 | measure | unaligned | passed | 1.929 | 0.773 | 0.842 | 0.704 | +19.7% |
| 256 x 256 x 64 | patient | matched | passed | 1.963 | 659.432 | 659.532 | 0.747 | +88200.7% |
| 256 x 256 x 64 | patient | unaligned | passed | 1.963 | 645.725 | 645.879 | 0.756 | +85352.6% |
| 256 x 256 x 64 | exhaustive | matched | passed | 1.796 | 641.617 | 641.693 | 0.771 | +83169.2% |
| 256 x 256 x 64 | exhaustive | unaligned | passed | 1.870 | 641.068 | 641.150 | 0.723 | +88558.6% |
| 1024 x 1024 x 30 | estimate | matched | passed | 18.474 | 9.538 | 9.645 | 9.510 | +1.4% |
| 1024 x 1024 x 30 | estimate | unaligned | passed | 18.618 | 9.968 | 10.052 | 9.923 | +1.3% |
| 1024 x 1024 x 30 | measure | matched | passed | 19.767 | 9.919 | 10.025 | 8.663 | +15.7% |
| 1024 x 1024 x 30 | measure | unaligned | passed | 22.079 | 10.190 | 10.302 | 8.537 | +20.7% |
| 1024 x 1024 x 30 | patient | matched | passed | 19.012 | 1582.436 | 1582.563 | 8.627 | +18244.8% |
| 1024 x 1024 x 30 | patient | unaligned | passed | 19.356 | 1668.838 | 1668.961 | 9.531 | +17411.5% |
| 1024 x 1024 x 30 | exhaustive | matched | passed | 19.652 | 1609.615 | 1609.731 | 8.924 | +17938.7% |
| 1024 x 1024 x 30 | exhaustive | unaligned | passed | 18.489 | 1572.572 | 1572.679 | 8.735 | +17904.3% |

## Readiness by workload

Speedups above 1 are faster than MATLAB. Positive bundled/native percentages mean bundled FFTW was slower.

| Size | Best bundled configuration | MATLAB forward (ms) | Bundled raw (ms) | Raw speedup | Bundled total (ms) | Total speedup | Bundled time relative to native | Destructive c2r speedup | Relative error | Half/full storage | Result |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| 256 x 256 x 64 | measure / unaligned / 18 threads | 1.929 | 0.773 | 2.495x | 0.842 | 2.292x | +19.7% | 2.173x | 5.02e-16 | 0.504 | READY |
| 1024 x 1024 x 30 | estimate / matched / 18 threads | 18.474 | 9.538 | 1.937x | 9.645 | 1.915x | +1.4% | 1.900x | 4.4e-16 | 0.501 | READY |

## Criterion status

| Size | Raw r2c >= 1.25x | Total MEX r2c >= 1.10x | Destructive c2r >= 0.95x | Error <= 1e-12 | Half spectrum | Zero copy |
|---|:---:|:---:|:---:|:---:|:---:|:---:|
| 256 x 256 x 64 | yes | yes | yes | yes | yes | yes |
| 1024 x 1024 x 30 | yes | yes | yes | yes | yes | yes |

## Timing boundaries

- Raw r2c time is the FFTW execute call only; no packing or layout conversion is required.
- Internal pipeline time includes MATLAB-buffer allocation, FFT execution, and buffer wrapping.
- Complete MEX time includes the entire MATLAB-to-MEX call. Input generation, planning, destructive-input refresh, metrics retrieval, and inverse normalization are excluded.
