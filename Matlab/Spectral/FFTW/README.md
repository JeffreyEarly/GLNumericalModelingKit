# FFTW transforms moved to FFTWTransforms

The reusable MATLAB FFTW backend formerly developed in this directory now
lives in the public
[`JeffreyEarly/fftw-transforms`](https://github.com/JeffreyEarly/fftw-transforms)
repository and is distributed through OceanKit as the `FFTWTransforms`
package. This directory is intentionally documentation-only so that adding
GLNumericalModelingKit to the MATLAB path does not expose duplicate transform
classes or MEX modules.

Use
[`FFTWTransforms 1.0.1`](https://github.com/JeffreyEarly/fftw-transforms/releases/tag/v1.0.1),
the verified corrective release. A package dependency of `^1.0.0` accepts this
release.

## Installation and local build

Register an OceanKit checkout and install the package by identity:

```matlab
mpmAddRepository("OceanKit","path/to/OceanKit")
mpminstall("FFTWTransforms")
```

The package distributes source only. It builds its MEX modules locally against
the FFTW library bundled with the active MATLAB installation; it does not
distribute FFTW libraries or compiled MEX files. Inspect support before asking
for a build:

```matlab
capabilities = FFTWBackend.capabilities();
if capabilities.build.isRequired && capabilities.build.isPossible
    capabilities = FFTWBackend.build();
end
```

The currently validated available backend is MATLAB R2026a on Apple silicon
(`maca64`). Supported package installations on other environments report
structured unavailability rather than attempting an unsupported compilation.
The verified OceanKit source snapshot was added in
[OceanKit commit `aa21cb7`](https://github.com/JeffreyEarly/OceanKit/commit/aa21cb7d455651a7d19092a66c71f6873b8f8a0d).

## Historical record

The retired GL implementation, experiments, tests, and canonical benchmark
artifacts remain available in Git history at
[GL commit `3c5fce2`](https://github.com/JeffreyEarly/GLNumericalModelingKit/tree/3c5fce2b2df9c892418676ef75ffbc5752216e55/Matlab/Spectral/FFTW).
The extracted package preserves the relevant commit history and its
[benchmark history](https://github.com/JeffreyEarly/fftw-transforms/tree/main/tools/benchmarks).

The GL work was developed and reviewed through these issues:

- [#37 — reproducible MATLAB baseline](https://github.com/JeffreyEarly/GLNumericalModelingKit/issues/37)
- [#38 — engine and half-spectrum layouts](https://github.com/JeffreyEarly/GLNumericalModelingKit/issues/38)
- [#39 — ownership, alignment, and copy costs](https://github.com/JeffreyEarly/GLNumericalModelingKit/issues/39)
- [#40 — backend support and readiness decision](https://github.com/JeffreyEarly/GLNumericalModelingKit/issues/40)
- [#41 — bundled versus native FFTW](https://github.com/JeffreyEarly/GLNumericalModelingKit/issues/41)
- [#42 — zero-copy half-spectrum backend](https://github.com/JeffreyEarly/GLNumericalModelingKit/issues/42)
- [#43 — DCT-I/DST-I backend](https://github.com/JeffreyEarly/GLNumericalModelingKit/issues/43)
- [#44 — release-aware build and capability detection](https://github.com/JeffreyEarly/GLNumericalModelingKit/issues/44)
- [#46 — retirement of the GL copy](https://github.com/JeffreyEarly/GLNumericalModelingKit/issues/46)

Issues [#32](https://github.com/JeffreyEarly/GLNumericalModelingKit/issues/32)
and [#34](https://github.com/JeffreyEarly/GLNumericalModelingKit/issues/34)
remain unchanged as earlier historical context.
