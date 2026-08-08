# Experimental FFTW transforms

This directory contains experimental MATLAB and MEX implementations for evaluating direct FFTW transforms. They are feasibility prototypes, not a supported GLNumericalModelingKit API or a WaveVortex integration commitment.

## Reproduce the issue #37 baseline

From the repository root, run:

```sh
matlab -batch "addpath('Matlab/Spectral/FFTW'); runFFTWFeasibilityBaseline"
```

The command rebuilds `fftw_dft2` against the active MATLAB installation's bundled FFTW library, runs the complete benchmark matrix, checks numerical and destructive-input behavior, and writes a new timestamped directory under `results/` containing `benchmark.json` and `summary.md`.

The benchmark records transform-call timing and memory facts without making the milestone's GO or NO-GO decision. MATLAB internal allocations and possible `complex(...)` or copy-on-write costs remain explicitly unresolved for issue #39.
