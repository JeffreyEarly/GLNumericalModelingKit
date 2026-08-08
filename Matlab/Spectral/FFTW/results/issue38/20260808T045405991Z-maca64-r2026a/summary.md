# FFT engine and half-spectrum layout benchmark

Status: **PASSED**

This issue #38 report shows criterion results but does not make the milestone's GO or NO-GO decision.

## Environment

| Field | Value |
|---|---|
| Run | 20260808T045405991Z-maca64-r2026a |
| MATLAB | 26.1.0.3312084 (R2026a) Update 4 |
| Architecture | maca64 |
| Operating system | Darwin 25.5.0 Darwin Kernel Version 25.5.0: Tue Jun  9 22:28:34 PDT 2026; root:xnu-12377.121.10~1/RELEASE_ARM64_T6050 arm64 |
| Processor | Apple M5 Max |
| Machine | MacBook Pro Mac17,6 |
| Physical memory | 48 GB |
| Threads | 18 |
| Planner limit | 10 seconds per plan |

## Engines

| Engine | Status | Version | Library | Reason |
|---|---|---|---|---|
| bundled-fftw | available | fftw-3.3.8 | /Applications/MATLAB_R2026a.app/bin/maca64/libmwfftw3.3.dylib |  |
| native-fftw | available | fftw-3.3.11 | /Users/jearly/Documents/OceanKitRepositories/GLNumericalModelingKit/Matlab/Spectral/FFTW/build/native-fftw-3.3.11/install/lib/libfftw3.3.dylib |  |
| accelerate-vdsp | available | Apple Accelerate/vDSP | /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libvDSP.dylib |  |

## FFTW new-array alignment validation

| Engine | Matched accepted | Mismatch rejected | Distinct offset classes observed |
|---|:---:|:---:|:---:|
| bundled-fftw | yes | yes | no |
| native-fftw | yes | yes | no |

## Fastest valid configuration by workload

Speed ratios above 1 are faster than the corresponding complete MATLAB transform.

| Size | Engine | Layout | Strategy | Planner | Threads | Alignment | Output | Raw r2c | Total MEX r2c | Destructive c2r | Relative error | Half/full |
|---|---|---|---|---|---:|---|---|---:|---:|---:|---:|---:|
| 128 x 128 x 128 | bundled-fftw | half-x | guru-rank2 | measure | 18 | matched | allocating | 2.316x | 1.565x | 1.042x | 3.48e-16 | 0.507812 |
| 256 x 128 x 40 | native-fftw | half-x | guru-rank2 | measure | 12 | unaligned | allocating | 2.395x | 1.721x | 1.314x | 4.4e-16 | 0.503906 |
| 256 x 256 x 64 | native-fftw | half-x | guru-rank2 | exhaustive | 18 | unaligned | allocating | 2.378x | 1.831x | 0.979x | 5.15e-16 | 0.503906 |
| 1024 x 1024 x 30 | native-fftw | half-x | guru-rank2 | patient | 18 | unaligned | allocating | 2.036x | 1.617x | 0.981x | 5.52e-16 | 0.500977 |

## Criterion status

The performance gate applies only to `256 x 256 x 64` and `1024 x 1024 x 30`.

| Size | Gate workload | Raw r2c >= 1.25x | Total MEX r2c >= 1.10x | Destructive c2r >= 0.95x | Error <= 1e-12 |
|---|:---:|:---:|:---:|:---:|:---:|
| 128 x 128 x 128 | no | yes | yes | yes | yes |
| 256 x 128 x 40 | no | yes | yes | yes | yes |
| 256 x 256 x 64 | yes | yes | yes | yes | yes |
| 1024 x 1024 x 30 | yes | yes | yes | yes | yes |

## Timing boundaries

- Kernel time contains only FFTW execute calls or vDSP transform calls.
- Raw pipeline time also contains required vDSP packing, unpacking, scaling, and canonical layout conversion.
- Total MEX time contains the complete MATLAB-to-MEX call using the existing allocation models. Detailed ownership attribution remains issue #39.
