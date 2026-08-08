# Experimental FFTW transforms

This directory contains experimental MATLAB and MEX implementations for evaluating direct FFTW transforms. They are feasibility prototypes, not a supported GLNumericalModelingKit API or a WaveVortex integration commitment.

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
