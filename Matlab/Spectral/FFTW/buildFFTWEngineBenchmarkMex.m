function build = buildFFTWEngineBenchmarkMex(options)
% Build the experimental issue #38 FFT engine benchmark gateways.
arguments
    options.shouldBuildNative (1,1) logical = false
    options.shouldBuildBundledAndVDSP (1,1) logical = true
    options.nativeBuildDirectory (1,1) string = fullfile(fileparts(mfilename('fullpath')),"build","native-fftw-3.3.11")
end

sourceDirectory = fileparts(mfilename('fullpath'));
includeArgument = "-I" + sourceDirectory;
bundledLibrary = fullfile(matlabroot,"bin",computer('arch'),"libmwfftw3.3.dylib");

build.bundled.module = "fftw_engine_benchmark_bundled";
build.bundled.library = string(bundledLibrary);
build.bundled.status = "existing";

build.vdsp.module = "vdsp_engine_benchmark";
build.vdsp.status = "existing";
build.vdsp.reason = "";
if options.shouldBuildBundledAndVDSP
    clear fftw_engine_benchmark_bundled
    mex('-R2018a','-outdir',sourceDirectory,'-output','fftw_engine_benchmark_bundled',includeArgument,fullfile(sourceDirectory,'fftw_engine_benchmark.cpp'),bundledLibrary);
    build.bundled.status = "built";
end
if options.shouldBuildBundledAndVDSP && ismac
    clear vdsp_engine_benchmark
    mex('-R2018a','-outdir',sourceDirectory,'-output','vdsp_engine_benchmark',fullfile(sourceDirectory,'vdsp_engine_benchmark.cpp'),'LDFLAGS=$LDFLAGS -framework Accelerate');
    build.vdsp.status = "built";
    build.vdsp.reason = "";
elseif ~ismac
    build.vdsp.status = "unavailable";
    build.vdsp.reason = "Accelerate/vDSP is available only on Apple platforms.";
end

build.native.module = "fftw_engine_benchmark_native";
build.native.status = "not-requested";
build.native.reason = "The bundled FFTW threshold result has not requested the native fallback.";
if options.shouldBuildNative
    native = buildNativeFFTW(options.nativeBuildDirectory);
    clear fftw_engine_benchmark_native
    nativeIncludeArgument = "-I" + native.includeDirectory;
    mex('-R2018a','-outdir',sourceDirectory,'-output','fftw_engine_benchmark_native',nativeIncludeArgument,fullfile(sourceDirectory,'fftw_engine_benchmark.cpp'),native.threadLibrary,native.baseLibrary);
    build.native = native;
    build.native.module = "fftw_engine_benchmark_native";
    build.native.status = "built";
    build.native.reason = "";
end
rehash
end

function native = buildNativeFFTW(buildDirectory)
version = "3.3.11";
archiveName = "fftw-" + version + ".tar.gz";
sourceURL = "https://fftw.org/pub/fftw/" + archiveName;
expectedSHA256 = "5630c24cdeb33b131612f7eb4b1a9934234754f9f388ff8617458d0be6f239a1";
configureFlags = "--enable-threads --disable-fortran --disable-openmp --enable-shared --disable-static";
deploymentTarget = "13.3";

archiveDirectory = fullfile(buildDirectory,"downloads");
sourceParent = fullfile(buildDirectory,"source");
installDirectory = fullfile(buildDirectory,"install");
archivePath = fullfile(archiveDirectory,archiveName);
sourceDirectory = fullfile(sourceParent,"fftw-" + version);

if ~isfolder(archiveDirectory), mkdir(archiveDirectory); end
if ~isfolder(sourceParent), mkdir(sourceParent); end
if ~isfile(archivePath)
    websave(archivePath,sourceURL);
end
actualSHA256 = binaryFileSHA256(archivePath);
if actualSHA256 ~= expectedSHA256
    error('FFTWEngineBenchmark:NativeChecksumMismatch','FFTW %s SHA-256 was %s; expected %s.',version,actualSHA256,expectedSHA256);
end
if ~isfolder(sourceDirectory)
    untar(archivePath,sourceParent);
end

baseLibrary = fullfile(installDirectory,"lib","libfftw3.3.dylib");
threadLibrary = fullfile(installDirectory,"lib","libfftw3_threads.3.dylib");
if ~isfile(baseLibrary) || ~isfile(threadLibrary)
    configureCommand = sprintf('cd "%s" && env MACOSX_DEPLOYMENT_TARGET="%s" CFLAGS="-O3 -mcpu=native -mmacosx-version-min=%s" LDFLAGS="-mmacosx-version-min=%s" ./configure --prefix="%s" %s',sourceDirectory,deploymentTarget,deploymentTarget,deploymentTarget,installDirectory,configureFlags);
    runBuildCommand(configureCommand,"configure FFTW");
    requestedJobs = max(1,maxNumCompThreads);
    runBuildCommand(sprintf('cd "%s" && make -j%d',sourceDirectory,requestedJobs),"compile FFTW");
    runBuildCommand(sprintf('cd "%s" && make install',sourceDirectory),"install FFTW");
end
if ~isfile(baseLibrary) || ~isfile(threadLibrary)
    error('FFTWEngineBenchmark:NativeLibraryMissing','The native FFTW build completed without producing the expected libraries.');
end

native.version = version;
native.sourceURL = sourceURL;
native.sourceSHA256 = actualSHA256;
native.configureFlags = configureFlags;
native.compilerFlags = "-O3 -mcpu=native -mmacosx-version-min=" + deploymentTarget + "; MACOSX_DEPLOYMENT_TARGET=" + deploymentTarget;
native.buildDirectory = buildDirectory;
native.installDirectory = string(installDirectory);
native.includeDirectory = string(fullfile(installDirectory,"include"));
native.baseLibrary = string(baseLibrary);
native.threadLibrary = string(threadLibrary);
end

function runBuildCommand(command,description)
[status,output] = system(command);
if status ~= 0
    error('FFTWEngineBenchmark:NativeBuildFailed','Unable to %s.\n%s',description,output);
end
end

function hash = binaryFileSHA256(path)
[fileId,message] = fopen(path,'r');
if fileId < 0
    error('FFTWEngineBenchmark:NativeArchiveOpenFailed','Unable to open %s: %s',path,message);
end
fileCleanup = onCleanup(@() fclose(fileId));
bytes = fread(fileId,Inf,'*uint8');
digest = java.security.MessageDigest.getInstance('SHA-256');
digest.update(bytes);
hashBytes = typecast(digest.digest(),'uint8');
hash = string(lower(reshape(dec2hex(hashBytes,2).',1,[])));
clear fileCleanup
end
