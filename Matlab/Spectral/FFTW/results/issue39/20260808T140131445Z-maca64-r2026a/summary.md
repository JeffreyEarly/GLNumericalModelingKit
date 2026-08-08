# FFTW MEX ownership, alignment, and copy benchmark

Status: **PASSED**

This issue #39 report characterizes ownership and copies without making the milestone's GO or NO-GO decision.

## Environment

| Field | Value |
|---|---|
| Run | 20260808T140131445Z-maca64-r2026a |
| MATLAB | 26.1.0.3312084 (R2026a) Update 4 |
| Architecture | maca64 |
| Operating system | Darwin 25.5.0 Darwin Kernel Version 25.5.0: Tue Jun  9 22:28:34 PDT 2026; root:xnu-12377.121.10~1/RELEASE_ARM64_T6050 arm64 arm |
| Processor | Apple M5 Max |
| Physical memory | 48.000 GB |

## Ownership timing

Times are medians in milliseconds. Boundary residual is complete MEX time minus the instrumented internal pipeline.

| Size | Operation | Total | Kernel | Allocate | Wrap | memcpy | Detach | Boundary | Detected copied MiB |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 128 x 128 x 128 | forward-factory-array | 0.4889 | 0.3907 | 0.0709 | 0.0000 | 0.0000 | 0.0000 | 0.0186 | 0.000 |
| 128 x 128 x 128 | forward-caller-direct | 0.4012 | 0.3752 | 0.0000 | 0.0000 | 0.0000 | 0.0002 | 0.0233 | 0.000 |
| 128 x 128 x 128 | forward-caller-complex | 0.4121 | 0.3945 | 0.0000 | 0.0000 | 0.0000 | 0.0001 | 0.0178 | 0.000 |
| 128 x 128 x 128 | forward-matlab-buffer | 0.3882 | 0.3700 | 0.0012 | 0.0022 | 0.0000 | 0.0000 | 0.0198 | 0.000 |
| 128 x 128 x 128 | forward-fftw-buffer | 0.4040 | 0.3765 | 0.0005 | 0.0034 | 0.0000 | 0.0000 | 0.0219 | 0.000 |
| 128 x 128 x 128 | inverse-preserving-allocating | 0.8032 | 0.4438 | 0.0713 | 0.0000 | 0.2558 | 0.0000 | 0.0248 | 0.000 |
| 128 x 128 x 128 | inverse-preserving-preallocated | 0.6827 | 0.4082 | 0.0000 | 0.0000 | 0.2522 | 0.0003 | 0.0238 | 0.000 |
| 128 x 128 x 128 | inverse-destructive-unique | 0.4687 | 0.4308 | 0.0000 | 0.0000 | 0.0000 | 0.0006 | 0.0372 | 0.000 |
| 128 x 128 x 128 | inverse-destructive-aliased | 0.7285 | 0.4435 | 0.0000 | 0.0000 | 0.0000 | 0.2397 | 0.0385 | 16.250 |
| 128 x 128 x 128 | inverse-destructive-complex | 0.4600 | 0.4316 | 0.0000 | 0.0000 | 0.0000 | 0.0005 | 0.0301 | 0.000 |
| 128 x 128 x 128 | chain-fftw-forward-destructive-inverse | 0.8588 | 0.7964 | 0.0009 | 0.0051 | 0.0000 | 0.0006 | 0.0500 | 0.000 |
| 256 x 128 x 40 | forward-factory-array | 0.3953 | 0.3054 | 0.0704 | 0.0000 | 0.0000 | 0.0000 | 0.0183 | 0.000 |
| 256 x 128 x 40 | forward-caller-direct | 0.3285 | 0.3100 | 0.0000 | 0.0000 | 0.0000 | 0.0002 | 0.0181 | 0.000 |
| 256 x 128 x 40 | forward-caller-complex | 0.2866 | 0.2696 | 0.0000 | 0.0000 | 0.0000 | 0.0001 | 0.0177 | 0.000 |
| 256 x 128 x 40 | forward-matlab-buffer | 0.3023 | 0.2826 | 0.0013 | 0.0015 | 0.0000 | 0.0000 | 0.0167 | 0.000 |
| 256 x 128 x 40 | forward-fftw-buffer | 0.3257 | 0.3041 | 0.0006 | 0.0035 | 0.0000 | 0.0000 | 0.0195 | 0.000 |
| 256 x 128 x 40 | inverse-preserving-allocating | 0.5827 | 0.3287 | 0.0721 | 0.0000 | 0.1573 | 0.0000 | 0.0208 | 0.000 |
| 256 x 128 x 40 | inverse-preserving-preallocated | 0.4867 | 0.3119 | 0.0000 | 0.0000 | 0.1544 | 0.0003 | 0.0208 | 0.000 |
| 256 x 128 x 40 | inverse-destructive-unique | 0.3460 | 0.3121 | 0.0000 | 0.0000 | 0.0000 | 0.0005 | 0.0333 | 0.000 |
| 256 x 128 x 40 | inverse-destructive-aliased | 0.4942 | 0.3120 | 0.0000 | 0.0000 | 0.0000 | 0.1470 | 0.0342 | 10.078 |
| 256 x 128 x 40 | inverse-destructive-complex | 0.3454 | 0.3198 | 0.0000 | 0.0000 | 0.0000 | 0.0005 | 0.0248 | 0.000 |
| 256 x 128 x 40 | chain-fftw-forward-destructive-inverse | 0.6897 | 0.6452 | 0.0011 | 0.0029 | 0.0000 | 0.0003 | 0.0398 | 0.000 |
| 256 x 256 x 64 | forward-factory-array | 1.3162 | 1.1213 | 0.1457 | 0.0000 | 0.0000 | 0.0000 | 0.0286 | 0.000 |
| 256 x 256 x 64 | forward-caller-direct | 1.1289 | 1.1013 | 0.0000 | 0.0000 | 0.0000 | 0.0003 | 0.0286 | 0.000 |
| 256 x 256 x 64 | forward-caller-complex | 1.1518 | 1.1206 | 0.0000 | 0.0000 | 0.0000 | 0.0002 | 0.0332 | 0.000 |
| 256 x 256 x 64 | forward-matlab-buffer | 1.1526 | 1.1132 | 0.0040 | 0.0040 | 0.0000 | 0.0000 | 0.0332 | 0.000 |
| 256 x 256 x 64 | forward-fftw-buffer | 1.1455 | 1.1094 | 0.0025 | 0.0064 | 0.0000 | 0.0000 | 0.0318 | 0.000 |
| 256 x 256 x 64 | inverse-preserving-allocating | 1.8362 | 1.1365 | 0.1452 | 0.0000 | 0.5401 | 0.0000 | 0.0307 | 0.000 |
| 256 x 256 x 64 | inverse-preserving-preallocated | 1.6187 | 1.0492 | 0.0000 | 0.0000 | 0.5234 | 0.0006 | 0.0368 | 0.000 |
| 256 x 256 x 64 | inverse-destructive-unique | 1.1998 | 1.1445 | 0.0000 | 0.0000 | 0.0000 | 0.0012 | 0.0583 | 0.000 |
| 256 x 256 x 64 | inverse-destructive-aliased | 1.7176 | 1.1354 | 0.0000 | 0.0000 | 0.0000 | 0.5168 | 0.0611 | 32.250 |
| 256 x 256 x 64 | inverse-destructive-complex | 1.1340 | 1.0875 | 0.0000 | 0.0000 | 0.0000 | 0.0008 | 0.0380 | 0.000 |
| 256 x 256 x 64 | chain-fftw-forward-destructive-inverse | 2.1578 | 2.0670 | 0.0026 | 0.0078 | 0.0000 | 0.0008 | 0.0701 | 0.000 |
| 1024 x 1024 x 30 | forward-factory-array | 9.3170 | 8.1429 | 1.1319 | 0.0000 | 0.0000 | 0.0000 | 0.0435 | 0.000 |
| 1024 x 1024 x 30 | forward-caller-direct | 9.1916 | 9.1445 | 0.0000 | 0.0000 | 0.0000 | 0.0004 | 0.0413 | 0.000 |
| 1024 x 1024 x 30 | forward-caller-complex | 8.6985 | 8.6413 | 0.0000 | 0.0000 | 0.0000 | 0.0004 | 0.0481 | 0.000 |
| 1024 x 1024 x 30 | forward-matlab-buffer | 8.0996 | 8.0179 | 0.0090 | 0.0076 | 0.0000 | 0.0000 | 0.0445 | 0.000 |
| 1024 x 1024 x 30 | forward-fftw-buffer | 9.1950 | 9.1295 | 0.0069 | 0.0142 | 0.0000 | 0.0000 | 0.0496 | 0.000 |
| 1024 x 1024 x 30 | inverse-preserving-allocating | 14.5266 | 9.6333 | 1.0556 | 0.0000 | 3.7606 | 0.0000 | 0.0465 | 0.000 |
| 1024 x 1024 x 30 | inverse-preserving-preallocated | 13.4209 | 9.6156 | 0.0000 | 0.0000 | 3.7395 | 0.0004 | 0.0481 | 0.000 |
| 1024 x 1024 x 30 | inverse-destructive-unique | 9.7137 | 9.6245 | 0.0000 | 0.0000 | 0.0000 | 0.0009 | 0.0853 | 0.000 |
| 1024 x 1024 x 30 | inverse-destructive-aliased | 12.0236 | 8.1505 | 0.0000 | 0.0000 | 0.0000 | 3.7798 | 0.0830 | 240.469 |
| 1024 x 1024 x 30 | inverse-destructive-complex | 8.6747 | 8.5956 | 0.0000 | 0.0000 | 0.0000 | 0.0007 | 0.0560 | 0.000 |
| 1024 x 1024 x 30 | chain-fftw-forward-destructive-inverse | 16.7405 | 16.6193 | 0.0073 | 0.0135 | 0.0000 | 0.0009 | 0.0976 | 0.000 |

### Allocation and copy ledger

| Size | Operation | Allocation owner | Timed allocated MiB | Explicit copied MiB | Detected COW MiB | Caller preallocated MiB | Persistent scratch MiB |
|---|---|---|---:|---:|---:|---:|---:|
| 128 x 128 x 128 | forward-factory-array | MATLAB Data API TypedArray | 16.250 | 0.000 | 0.000 | 0.000 | 16.250 |
| 128 x 128 x 128 | forward-caller-direct | caller; MATLAB copy-on-write only if detached | 0.000 | 0.000 | 0.000 | 16.250 | 16.250 |
| 128 x 128 x 128 | forward-caller-complex | caller; MATLAB copy-on-write only if detached | 0.000 | 0.000 | 0.000 | 16.250 | 16.250 |
| 128 x 128 x 128 | forward-matlab-buffer | MATLAB Data API buffer | 16.250 | 0.000 | 0.000 | 0.000 | 16.250 |
| 128 x 128 x 128 | forward-fftw-buffer | FFTW buffer with custom deleter | 16.250 | 0.000 | 0.000 | 0.000 | 16.250 |
| 128 x 128 x 128 | inverse-preserving-allocating | MATLAB Data API real output plus persistent MEX scratch | 16.000 | 16.250 | 0.000 | 0.000 | 16.250 |
| 128 x 128 x 128 | inverse-preserving-preallocated | caller real output plus persistent MEX scratch | 0.000 | 16.250 | 0.000 | 16.000 | 16.250 |
| 128 x 128 x 128 | inverse-destructive-unique | caller spectrum and real output; MATLAB copy-on-write if detached | 0.000 | 0.000 | 0.000 | 32.250 | 16.250 |
| 128 x 128 x 128 | inverse-destructive-aliased | caller spectrum and real output; MATLAB copy-on-write if detached | 16.250 | 0.000 | 16.250 | 32.250 | 16.250 |
| 128 x 128 x 128 | inverse-destructive-complex | caller spectrum and real output; MATLAB copy-on-write if detached | 0.000 | 0.000 | 0.000 | 32.250 | 16.250 |
| 128 x 128 x 128 | chain-fftw-forward-destructive-inverse | FFTW spectrum with custom deleter plus caller real output | 16.250 | 0.000 | 0.000 | 16.000 | 16.250 |
| 256 x 128 x 40 | forward-factory-array | MATLAB Data API TypedArray | 10.078 | 0.000 | 0.000 | 0.000 | 10.078 |
| 256 x 128 x 40 | forward-caller-direct | caller; MATLAB copy-on-write only if detached | 0.000 | 0.000 | 0.000 | 10.078 | 10.078 |
| 256 x 128 x 40 | forward-caller-complex | caller; MATLAB copy-on-write only if detached | 0.000 | 0.000 | 0.000 | 10.078 | 10.078 |
| 256 x 128 x 40 | forward-matlab-buffer | MATLAB Data API buffer | 10.078 | 0.000 | 0.000 | 0.000 | 10.078 |
| 256 x 128 x 40 | forward-fftw-buffer | FFTW buffer with custom deleter | 10.078 | 0.000 | 0.000 | 0.000 | 10.078 |
| 256 x 128 x 40 | inverse-preserving-allocating | MATLAB Data API real output plus persistent MEX scratch | 10.000 | 10.078 | 0.000 | 0.000 | 10.078 |
| 256 x 128 x 40 | inverse-preserving-preallocated | caller real output plus persistent MEX scratch | 0.000 | 10.078 | 0.000 | 10.000 | 10.078 |
| 256 x 128 x 40 | inverse-destructive-unique | caller spectrum and real output; MATLAB copy-on-write if detached | 0.000 | 0.000 | 0.000 | 20.078 | 10.078 |
| 256 x 128 x 40 | inverse-destructive-aliased | caller spectrum and real output; MATLAB copy-on-write if detached | 10.078 | 0.000 | 10.078 | 20.078 | 10.078 |
| 256 x 128 x 40 | inverse-destructive-complex | caller spectrum and real output; MATLAB copy-on-write if detached | 0.000 | 0.000 | 0.000 | 20.078 | 10.078 |
| 256 x 128 x 40 | chain-fftw-forward-destructive-inverse | FFTW spectrum with custom deleter plus caller real output | 10.078 | 0.000 | 0.000 | 10.000 | 10.078 |
| 256 x 256 x 64 | forward-factory-array | MATLAB Data API TypedArray | 32.250 | 0.000 | 0.000 | 0.000 | 32.250 |
| 256 x 256 x 64 | forward-caller-direct | caller; MATLAB copy-on-write only if detached | 0.000 | 0.000 | 0.000 | 32.250 | 32.250 |
| 256 x 256 x 64 | forward-caller-complex | caller; MATLAB copy-on-write only if detached | 0.000 | 0.000 | 0.000 | 32.250 | 32.250 |
| 256 x 256 x 64 | forward-matlab-buffer | MATLAB Data API buffer | 32.250 | 0.000 | 0.000 | 0.000 | 32.250 |
| 256 x 256 x 64 | forward-fftw-buffer | FFTW buffer with custom deleter | 32.250 | 0.000 | 0.000 | 0.000 | 32.250 |
| 256 x 256 x 64 | inverse-preserving-allocating | MATLAB Data API real output plus persistent MEX scratch | 32.000 | 32.250 | 0.000 | 0.000 | 32.250 |
| 256 x 256 x 64 | inverse-preserving-preallocated | caller real output plus persistent MEX scratch | 0.000 | 32.250 | 0.000 | 32.000 | 32.250 |
| 256 x 256 x 64 | inverse-destructive-unique | caller spectrum and real output; MATLAB copy-on-write if detached | 0.000 | 0.000 | 0.000 | 64.250 | 32.250 |
| 256 x 256 x 64 | inverse-destructive-aliased | caller spectrum and real output; MATLAB copy-on-write if detached | 32.250 | 0.000 | 32.250 | 64.250 | 32.250 |
| 256 x 256 x 64 | inverse-destructive-complex | caller spectrum and real output; MATLAB copy-on-write if detached | 0.000 | 0.000 | 0.000 | 64.250 | 32.250 |
| 256 x 256 x 64 | chain-fftw-forward-destructive-inverse | FFTW spectrum with custom deleter plus caller real output | 32.250 | 0.000 | 0.000 | 32.000 | 32.250 |
| 1024 x 1024 x 30 | forward-factory-array | MATLAB Data API TypedArray | 240.469 | 0.000 | 0.000 | 0.000 | 240.469 |
| 1024 x 1024 x 30 | forward-caller-direct | caller; MATLAB copy-on-write only if detached | 0.000 | 0.000 | 0.000 | 240.469 | 240.469 |
| 1024 x 1024 x 30 | forward-caller-complex | caller; MATLAB copy-on-write only if detached | 0.000 | 0.000 | 0.000 | 240.469 | 240.469 |
| 1024 x 1024 x 30 | forward-matlab-buffer | MATLAB Data API buffer | 240.469 | 0.000 | 0.000 | 0.000 | 240.469 |
| 1024 x 1024 x 30 | forward-fftw-buffer | FFTW buffer with custom deleter | 240.469 | 0.000 | 0.000 | 0.000 | 240.469 |
| 1024 x 1024 x 30 | inverse-preserving-allocating | MATLAB Data API real output plus persistent MEX scratch | 240.000 | 240.469 | 0.000 | 0.000 | 240.469 |
| 1024 x 1024 x 30 | inverse-preserving-preallocated | caller real output plus persistent MEX scratch | 0.000 | 240.469 | 0.000 | 240.000 | 240.469 |
| 1024 x 1024 x 30 | inverse-destructive-unique | caller spectrum and real output; MATLAB copy-on-write if detached | 0.000 | 0.000 | 0.000 | 480.469 | 240.469 |
| 1024 x 1024 x 30 | inverse-destructive-aliased | caller spectrum and real output; MATLAB copy-on-write if detached | 240.469 | 0.000 | 240.469 | 480.469 | 240.469 |
| 1024 x 1024 x 30 | inverse-destructive-complex | caller spectrum and real output; MATLAB copy-on-write if detached | 0.000 | 0.000 | 0.000 | 480.469 | 240.469 |
| 1024 x 1024 x 30 | chain-fftw-forward-destructive-inverse | FFTW spectrum with custom deleter plus caller real output | 240.469 | 0.000 | 0.000 | 240.000 | 240.469 |

## Zero-copy and lifetime evidence

| Size | Operation | Input detached | Output detached | Wrapper allocated | Return pointer preserved | Explicit copied MiB | Error |
|---|---|:---:|:---:|:---:|:---:|---:|---:|
| 128 x 128 x 128 | forward-factory-array | no | no | no | yes | 0.000 | 4.54e-16 |
| 128 x 128 x 128 | forward-caller-direct | no | no | no | yes | 0.000 | 4.54e-16 |
| 128 x 128 x 128 | forward-caller-complex | no | no | no | yes | 0.000 | 4.54e-16 |
| 128 x 128 x 128 | forward-matlab-buffer | no | no | no | yes | 0.000 | 4.54e-16 |
| 128 x 128 x 128 | forward-fftw-buffer | no | no | no | yes | 0.000 | 4.54e-16 |
| 128 x 128 x 128 | inverse-preserving-allocating | no | no | no | yes | 16.250 | 3.74e-16 |
| 128 x 128 x 128 | inverse-preserving-preallocated | no | no | no | yes | 16.250 | 3.74e-16 |
| 128 x 128 x 128 | inverse-destructive-unique | no | no | no | yes | 0.000 | 3.74e-16 |
| 128 x 128 x 128 | inverse-destructive-aliased | yes | no | no | yes | 0.000 | 3.74e-16 |
| 128 x 128 x 128 | inverse-destructive-complex | no | no | no | yes | 0.000 | 3.74e-16 |
| 128 x 128 x 128 | chain-fftw-forward-destructive-inverse | no | no | no | yes | 0.000 | 3.74e-16 |
| 256 x 128 x 40 | forward-factory-array | no | no | no | yes | 0.000 | 3.63e-16 |
| 256 x 128 x 40 | forward-caller-direct | no | no | no | yes | 0.000 | 3.63e-16 |
| 256 x 128 x 40 | forward-caller-complex | no | no | no | yes | 0.000 | 3.63e-16 |
| 256 x 128 x 40 | forward-matlab-buffer | no | no | no | yes | 0.000 | 3.63e-16 |
| 256 x 128 x 40 | forward-fftw-buffer | no | no | no | yes | 0.000 | 3.63e-16 |
| 256 x 128 x 40 | inverse-preserving-allocating | no | no | no | yes | 10.078 | 4.2e-16 |
| 256 x 128 x 40 | inverse-preserving-preallocated | no | no | no | yes | 10.078 | 4.2e-16 |
| 256 x 128 x 40 | inverse-destructive-unique | no | no | no | yes | 0.000 | 4.2e-16 |
| 256 x 128 x 40 | inverse-destructive-aliased | yes | no | no | yes | 0.000 | 4.2e-16 |
| 256 x 128 x 40 | inverse-destructive-complex | no | no | no | yes | 0.000 | 4.2e-16 |
| 256 x 128 x 40 | chain-fftw-forward-destructive-inverse | no | no | no | yes | 0.000 | 4.2e-16 |
| 256 x 256 x 64 | forward-factory-array | no | no | no | yes | 0.000 | 4.12e-16 |
| 256 x 256 x 64 | forward-caller-direct | no | no | no | yes | 0.000 | 4.12e-16 |
| 256 x 256 x 64 | forward-caller-complex | no | no | no | yes | 0.000 | 4.12e-16 |
| 256 x 256 x 64 | forward-matlab-buffer | no | no | no | yes | 0.000 | 4.12e-16 |
| 256 x 256 x 64 | forward-fftw-buffer | no | no | no | yes | 0.000 | 4.12e-16 |
| 256 x 256 x 64 | inverse-preserving-allocating | no | no | no | yes | 32.250 | 5.09e-16 |
| 256 x 256 x 64 | inverse-preserving-preallocated | no | no | no | yes | 32.250 | 5.09e-16 |
| 256 x 256 x 64 | inverse-destructive-unique | no | no | no | yes | 0.000 | 5.09e-16 |
| 256 x 256 x 64 | inverse-destructive-aliased | yes | no | no | yes | 0.000 | 5.09e-16 |
| 256 x 256 x 64 | inverse-destructive-complex | no | no | no | yes | 0.000 | 5.09e-16 |
| 256 x 256 x 64 | chain-fftw-forward-destructive-inverse | no | no | no | yes | 0.000 | 4.67e-16 |
| 1024 x 1024 x 30 | forward-factory-array | no | no | no | yes | 0.000 | 4.41e-16 |
| 1024 x 1024 x 30 | forward-caller-direct | no | no | no | yes | 0.000 | 4.41e-16 |
| 1024 x 1024 x 30 | forward-caller-complex | no | no | no | yes | 0.000 | 4.41e-16 |
| 1024 x 1024 x 30 | forward-matlab-buffer | no | no | no | yes | 0.000 | 4.41e-16 |
| 1024 x 1024 x 30 | forward-fftw-buffer | no | no | no | yes | 0.000 | 4.41e-16 |
| 1024 x 1024 x 30 | inverse-preserving-allocating | no | no | no | yes | 240.469 | 4.89e-16 |
| 1024 x 1024 x 30 | inverse-preserving-preallocated | no | no | no | yes | 240.469 | 4.89e-16 |
| 1024 x 1024 x 30 | inverse-destructive-unique | no | no | no | yes | 0.000 | 4.89e-16 |
| 1024 x 1024 x 30 | inverse-destructive-aliased | yes | no | no | yes | 0.000 | 4.89e-16 |
| 1024 x 1024 x 30 | inverse-destructive-complex | no | no | no | yes | 0.000 | 4.89e-16 |
| 1024 x 1024 x 30 | chain-fftw-forward-destructive-inverse | no | no | no | yes | 0.000 | 5.26e-16 |

| Size | FFTW buffers created | Freed | Outstanding | Balanced |
|---|---:|---:|---:|:---:|
| 128 x 128 x 128 | 18 | 18 | 0 | yes |
| 256 x 128 x 40 | 18 | 18 | 0 | yes |
| 256 x 256 x 64 | 18 | 18 | 0 | yes |
| 1024 x 1024 x 30 | 10 | 10 | 0 | yes |

## Persistent storage

| Size | MATLAB full spectrum MiB | Half spectrum MiB | Legacy standalone MiB | WaveVortex MATLAB MiB | Current FFTW backend MiB | Lean preserving scratch MiB | Estimated plan MiB |
|---|---:|---:|---:|---:|---:|---:|---:|
| 128 x 128 x 128 | 32.000 | 16.250 | 32.500 | 32.000 | 194.250 | 16.250 | 0.102 |
| 256 x 128 x 40 | 20.000 | 10.078 | 20.156 | 20.000 | 121.172 | 10.078 | 5.993 |
| 256 x 256 x 64 | 64.000 | 32.250 | 64.500 | 64.000 | 386.250 | 32.250 | 0.100 |
| 1024 x 1024 x 30 | 480.000 | 240.469 | 480.938 | 480.000 | 2884.219 | 240.469 | 0.105 |

Plan-owned values are allocator-delta estimates; application-owned buffers and copies are exact.

## Recommendation

- Forward ownership: `forward-matlab-buffer`.
- Forward contract: Allocate a MATLAB-managed buffer inside MEX and return it through createArrayFromBuffer.
- Selection basis: raw complete MEX medians.
- A common raw 5% candidate existed: **yes**.
- Destructive inverse demonstrated without a spectrum copy: **yes**.
- Contract: Move a uniquely owned half-spectrum directly, return and reassign the destroyed spectrum, and omit complex(...).
- Preserving inverse: Input-preserving c2r retains exactly one explicit half-spectrum memcpy into persistent scratch.