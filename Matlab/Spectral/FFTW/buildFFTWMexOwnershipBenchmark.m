function build = buildFFTWMexOwnershipBenchmark(options)
% Build the experimental issue #39 FFTW ownership benchmark gateways.
arguments
    options.shouldBuildNative (1,1) logical = true
    options.nativeBuildDirectory (1,1) string = fullfile(fileparts(mfilename('fullpath')),"build","native-fftw-3.3.11")
end

sourceDirectory = fileparts(mfilename('fullpath'));
sourcePath = fullfile(sourceDirectory,"fftw_ownership_benchmark.cpp");
bundledLibrary = fullfile(matlabroot,"bin",computer('arch'),"libmwfftw3.3.dylib");
includeArgument = "-I" + sourceDirectory;

clear fftw_ownership_benchmark_bundled
mex('-R2018a','-outdir',sourceDirectory,'-output','fftw_ownership_benchmark_bundled',includeArgument,sourcePath,bundledLibrary);
build.bundled.module = "fftw_ownership_benchmark_bundled";
build.bundled.library = string(bundledLibrary);
build.bundled.status = "built";

build.native.module = "fftw_ownership_benchmark_native";
build.native.status = "not-requested";
build.native.reason = "Native FFTW was not requested.";
if options.shouldBuildNative
    engineBuild = buildFFTWEngineBenchmarkMex(shouldBuildNative=true,shouldBuildBundledAndVDSP=false,nativeBuildDirectory=options.nativeBuildDirectory);
    native = engineBuild.native;
    clear fftw_ownership_benchmark_native
    nativeIncludeArgument = "-I" + native.includeDirectory;
    mex('-R2018a','-outdir',sourceDirectory,'-output','fftw_ownership_benchmark_native',nativeIncludeArgument,sourcePath,native.threadLibrary,native.baseLibrary);
    build.native = native;
    build.native.module = "fftw_ownership_benchmark_native";
    build.native.status = "built";
    build.native.reason = "";
end
rehash
end
