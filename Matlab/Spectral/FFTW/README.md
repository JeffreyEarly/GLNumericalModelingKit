# FFTW transforms

This directory contains the bundled-FFTW real-to-complex backend and the feasibility benchmarks that established its performance and ownership contracts. The benchmark gateways remain experimental and are not a WaveVortex integration commitment.

## Discover and build the backend

`FFTWBackend` is the release-aware entry point for normal use. A capability query never builds or warns, and expected availability failures are returned as structured status:

```matlab
capabilities = FFTWBackend.capabilities();
if capabilities.build.isRequired && capabilities.build.isPossible
    capabilities = FFTWBackend.build();
end
```

The initial provider is `matlab-bundled`, supported on MATLAB R2026a updates for macOS `maca64`. The builder resolves the active installation's FFTW library, compiles both production gateways in a temporary directory, validates their loaded library, alignment, and numerical behavior, and installs both only after the staged modules pass. A failed replacement restores the previous ignored MEX files.

The result has top-level `status`, `isAvailable`, and `isComplete` fields. Inspect `features.r2c`, `features.c2r`, `features.dct1`, and `features.dst1` individually; a partial backend can have only some usable features. Each module and feature contains a structured `reason` with a stable code, stage, module, identifier, and actionable message. The `eligibility` field contains the issue #41 horizontal readiness contract and issue #43's bounded DCT-I/DST-I table. Real-to-real eligibility applies only to an exact tested `Nz` and an inclusive recorded batch interval.

Unsupported releases, architectures, missing compilers or libraries, locked modules, compilation failures, library mismatches, and failed self-tests are returned rather than thrown. No FFTW source, FFTW library, or precompiled MEX binary is distributed.

## Build and use the half-spectrum backend

Build the production gateway against the active MATLAB installation's bundled FFTW library from any working directory:

```matlab
addpath('/path/to/GLNumericalModelingKit/Matlab/Spectral/FFTW')
FFTWBackend.build()
```

`RealToComplexTransform.makeMexFiles` remains available as a compatibility build entry point.

`RealToComplexTransform` accepts a nonempty ordered list of distinct transform dimensions. The final entry is the compressed dimension: `[2 1]` produces half-x storage, `[1 2]` produces half-y storage, and `[3 1]` transforms z before storing half-x. If the real length of compressed dimension `d` is `N`, the complex output length in `d` is `floor(N/2)+1`; all other lengths are unchanged.

```matlab
transform = RealToComplexTransform([128 128 32],dims=[2 1]);
spectrum = transform.transformForward(realArray);
roundTrip = transform.scaleFactor*transform.transformBack(spectrum);
```

FFTW's inverse is deliberately left unnormalized. `scaleFactor` is `1/prod(transform.realSize(transform.transformDimensions))`, allowing a caller to fuse normalization with later spectral operations.

Allocating forward transforms return a MATLAB-managed half-spectrum buffer without a spectrum copy. Preserving inverse methods perform exactly one half-spectrum copy because multidimensional FFTW c2r execution destroys its input. `transformBackIntoArrayDestructive` avoids that explicit copy for uniquely owned arrays and returns both the destroyed spectrum and real output. Always reassign caller-preallocated and destructive outputs. If either destructive input has an alias, MATLAB may detach it through copy-on-write; the aliases remain unchanged, but that call is outside the zero-copy guarantee.

The default `alignmentMode="unaligned"` accepts arbitrary MATLAB array alignment. `alignmentMode="matched"` is an explicit expert option: every execution must match the input and output alignment classes used to construct the plan, and incompatible arrays are rejected.

Run the production test suite with:

```sh
matlab -batch "addpath('Matlab/Spectral/FFTW'); results=runtests('Matlab/Spectral/FFTW/UnitTests/TestRealToComplexTransform.m'); assertSuccess(results)"
```

## Build and use the real-to-real backend

RealToRealTransform is the production DCT-I/DST-I interface and follows the same constructor, ownership, alignment, and build conventions as RealToComplexTransform:

~~~matlab
FFTWBackend.build()
transform = RealToRealTransform([65 33024],dims=1,transform="cosine",dataType="complex");
coefficients = transform.transformForward(values);
roundTrip = transform.transformBack(coefficients);
~~~

Cosine transforms retain the physical shape. Sine transforms omit the two physical endpoints from spectralSize; forward execution ignores their values and inverse execution restores exact zeros. Both directions are normalized to match the GL transform matrices, so scaleFactor is one. Complex arrays remain interleaved and are transformed as batched real and imaginary components without MATLAB packing arrays.

Caller-preallocated outputs must be reassigned. The default unaligned plans accept arbitrary MATLAB arrays; matched alignment remains an explicit expert option.

`RealToRealTransform.makeMexFiles` remains available as a compatibility build entry point.

Run the issue #43 benchmark with:

~~~sh
matlab -batch "addpath('Matlab/Spectral'); addpath('Matlab/Spectral/FFTW'); runFFTWR2RBenchmark"
~~~

The benchmark compares complete allocating calls with the dense matrices and FFT-extension implementations at the validated vertical sizes and half-x batch anchors. It writes raw timings, winners, and bounded eligibility intervals beneath results/issue43/.

Run the production and benchmark tests with:

~~~sh
matlab -batch "addpath('Matlab/Spectral'); addpath('Matlab/Spectral/FFTW'); results=runtests({'Matlab/Spectral/FFTW/UnitTests/TestRealToRealTransform.m','Matlab/Spectral/FFTW/UnitTests/TestFFTWR2RBenchmark.m'}); assertSuccess(results)"
~~~

## Reproduce the issue #37 baseline

From the repository root, run:

```sh
matlab -batch "addpath('Matlab/Spectral/FFTW'); runFFTWFeasibilityBaseline"
```

The command rebuilds `fftw_dft2` against the active MATLAB installation's bundled FFTW library, runs the complete benchmark matrix, checks numerical and destructive-input behavior, and writes a new timestamped directory under `results/` containing `benchmark.json` and `summary.md`.

The benchmark records transform-call timing and memory facts without making the milestone's GO or NO-GO decision. MATLAB internal allocations and possible `complex(...)` or copy-on-write costs remain explicitly unresolved for issue #39.

## Reproduce the issue #38 engine and layout benchmark

From the repository root, run:

```sh
matlab -batch "addpath('Matlab/Spectral/FFTW'); runFFTWEngineLayoutBenchmark"
```

The issue #38 benchmark first screens MATLAB's bundled FFTW across half-x and half-y layouts, rank-2 and staged transforms, all four FFTW planners, matched and unaligned execution, and thread counts up to the machine maximum. It then remeasures the fastest valid configuration for each workload with MATLAB operations interleaved in round-robin order.

If bundled FFTW misses the forward thresholds on either gate workload, the same command builds FFTW 3.3.11 from its pinned official source in the ignored `build/native-fftw-3.3.11` directory and also screens Accelerate/vDSP. Downloaded FFTW source and GPL-licensed native binaries are not committed. Each run writes `engine-layout-benchmark.json` and `summary.md` beneath a timestamped directory in `results/issue38/`.

The internal `kernelSeconds` measurement includes vendor transform calls only. `rawPipelineSeconds` additionally includes packing, unpacking, scaling, and canonical layout conversion required by the engine. `totalMexSeconds` covers the complete MATLAB-to-MEX call. Results show the issue #38 criteria but leave the milestone's GO or NO-GO decision to the decision issue.

Run the issue #38 smoke suite with:

```sh
matlab -batch "addpath('Matlab/Spectral/FFTW'); results=runtests('Matlab/Spectral/FFTW/UnitTests/TestFFTWEngineLayoutBenchmark.m'); assertSuccess(results)"
```

## Reproduce the issue #39 ownership and copy benchmark

From the repository root, run:

```sh
matlab -batch "addpath('Matlab/Spectral/FFTW'); runFFTWMexOwnershipBenchmark"
```

The issue #39 benchmark imports the committed issue #38 winner for each workload and separates FFT execution, output allocation, buffer wrapping, explicit `memcpy`, mutable-array detachment, and MATLAB/MEX boundary time. It compares ordinary MATLAB-owned arrays, caller-preallocated arrays, MATLAB Data API buffers, and FFTW-owned buffers returned with a custom `fftw_free` deleter. Pointer tokens and lifetime counters verify whether each path actually crosses the MEX boundary without a spectrum copy.

Each run writes `ownership-benchmark.json` and `summary.md` beneath a timestamped directory in `results/issue39/`. The report recommends an ownership contract but leaves the milestone's GO or NO-GO decision to issue #40.

Run the issue #39 smoke suite with:

```sh
matlab -batch "addpath('Matlab/Spectral/FFTW'); results=runtests('Matlab/Spectral/FFTW/UnitTests/TestFFTWMexOwnershipBenchmark.m'); assertSuccess(results)"
```
