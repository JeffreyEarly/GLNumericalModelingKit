function build = buildFFTWBundledNativeComparisonMex(options)
% Build the private issue #41 bundled/native comparison gateways.
arguments
    options.nativeBuildDirectory (1,1) string = fullfile(fileparts(mfilename('fullpath')),"build","native-fftw-3.3.11")
end

sourceDirectory = fileparts(mfilename('fullpath'));
build.ownership = buildFFTWMexOwnershipBenchmark(shouldBuildNative=true,nativeBuildDirectory=options.nativeBuildDirectory);

bundledLibrary = fullfile(matlabroot,"bin",computer('arch'),"libmwfftw3.3.dylib");
includeArgument = "-I" + sourceDirectory;
clear fftw_engine_benchmark_bundled
mex('-R2018a','-outdir',sourceDirectory,'-output','fftw_engine_benchmark_bundled',includeArgument,fullfile(sourceDirectory,'fftw_engine_benchmark.cpp'),bundledLibrary);

build.replay.bundled.module = "fftw_engine_benchmark_bundled";
build.replay.bundled.library = string(bundledLibrary);
build.replay.bundled.status = "built";
build.replay.native.module = "fftw_engine_benchmark_native";
build.replay.native.status = build.ownership.native.status;
build.replay.native.library = build.ownership.native.baseLibrary;
rehash
end
