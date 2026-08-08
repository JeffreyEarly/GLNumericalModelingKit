# FFTW r2c feasibility baseline

Status: **PASSED**

This report is a descriptive Stage 1 baseline. It does not make a GO or NO-GO decision.

## Environment

| Field | Value |
|---|---|
| Run | 20260808T040109179Z-maca64-r2026a |
| Generated (UTC) | 2026-08-08T04:01:09.231Z |
| MATLAB | 26.1.0.3312084 (R2026a) Update 4 |
| Architecture | maca64 |
| Operating system | Darwin 25.5.0 Darwin Kernel Version 25.5.0: Tue Jun  9 22:28:34 PDT 2026; root:xnu-12377.121.10~1/RELEASE_ARM64_T6050 arm64 |
| Processor | Apple M5 Max |
| Machine | MacBook Pro Mac17,6 |
| Physical memory | 48 GB |
| FFTW | fftw-3.3.8 via /Applications/MATLAB_R2026a.app/bin/maca64/libmwfftw3.3.dylib |
| Planner | measure |
| Threads | 18 |
| Warmups | 2 |
| Samples | 7 (3 for largest) |

## Median transform-call time

MEX cells show seconds followed by speed relative to the corresponding MATLAB transform. Ratios above 1 are faster than MATLAB.

| Size | MATLAB forward | MEX r2c alloc | MEX r2c prealloc | MATLAB inverse | MEX c2r alloc | MEX c2r prealloc | MEX c2r destructive |
|---|---:|---:|---:|---:|---:|---:|---:|
| 128 x 128 x 128 | 0.00070975 | 0.000779 (0.911x) | 0.000670292 (1.059x) | 0.000872083 | 0.000947792 (0.920x) | 0.000881333 (0.990x) | 0.000662958 (1.315x) |
| 256 x 128 x 40 | 0.00057475 | 0.000642792 (0.894x) | 0.000580875 (0.989x) | 0.000638625 | 0.000912167 (0.700x) | 0.000828583 (0.771x) | 0.000703417 (0.908x) |
| 256 x 256 x 64 | 0.00163958 | 0.00236263 (0.694x) | 0.00221125 (0.741x) | 0.00181125 | 0.00344183 (0.526x) | 0.00330617 (0.548x) | 0.00282142 (0.642x) |
| 1024 x 1024 x 30 | 0.0172285 | 0.030053 (0.573x) | 0.0309603 (0.556x) | 0.0172383 | 0.037978 (0.454x) | 0.0380425 (0.453x) | 0.0332568 (0.518x) |

## Correctness

The pass threshold is a relative infinity error of `1e-12`.

| Size | Maximum absolute error | Maximum relative error | Preserving input unchanged | Destructive input changed | Pass |
|---|---:|---:|:---:|:---:|:---:|
| 128 x 128 x 128 | 1.73502e-13 | 4.15485e-16 | yes | yes | yes |
| 256 x 128 x 40 | 2.84217e-13 | 4.19541e-16 | yes | yes | yes |
| 256 x 256 x 64 | 3.81317e-13 | 5.09381e-16 | yes | yes | yes |
| 1024 x 1024 x 30 | 1.94268e-12 | 4.53535e-16 | yes | yes | yes |

## Array storage

| Size | Real array (MiB) | MATLAB full spectrum (MiB) | MEX half spectrum (MiB) | Half/full ratio | Persistent transform buffers (MiB) |
|---|---:|---:|---:|---:|---:|
| 128 x 128 x 128 | 16.000 | 32.000 | 16.250 | 0.507812 | 32.500 |
| 256 x 128 x 40 | 10.000 | 20.000 | 10.156 | 0.507812 | 20.312 |
| 256 x 256 x 64 | 32.000 | 64.000 | 32.250 | 0.503906 | 64.500 |
| 1024 x 1024 x 30 | 240.000 | 480.000 | 240.469 | 0.500977 | 480.938 |

## Allocation model

Counts are source-grounded minimums. MATLAB internal FFT work buffers and possible `complex(...)` or copy-on-write allocations remain unresolved for issue #39.

| Operation | Minimum known timed allocation | Minimum known timed copy | Caller preallocation |
|---|---|---|---|
| MATLAB forward | Two full-spectrum arrays | Unknown internal work | None |
| MATLAB inverse | One full-spectrum and one real array | Unknown internal work | None |
| MEX r2c allocating | One half-spectrum output | None known | None |
| MEX r2c preallocated | None known | Copy-on-write unresolved | Half-spectrum output |
| MEX c2r allocating preserving | One real output | One half-spectrum `memcpy` | None |
| MEX c2r preallocated preserving | None known | One half-spectrum `memcpy` | Real output |
| MEX c2r preallocated destructive | None known | Copy-on-write unresolved | Half-spectrum input and real output |
