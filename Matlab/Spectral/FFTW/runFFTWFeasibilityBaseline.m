function result = runFFTWFeasibilityBaseline(options)
% Build and benchmark the experimental FFTW real-to-complex MEX backend.
%
% The default run exercises representative WaveVortex array sizes, writes
% machine-readable JSON and a Markdown summary, and throws if compilation or
% numerical validation fails. Transform planning, allocation setup, warmups,
% and normalization are excluded from the reported timings.
arguments
    options.outputDirectory (1,1) string = fullfile(fileparts(mfilename('fullpath')),"results")
    options.sizes (:,3) double = [128 128 128; 256 128 40; 256 256 64; 1024 1024 30]
    options.planner (1,1) string = "measure"
    options.nThreads (1,1) double = NaN
    options.nWarmups (1,1) double = 2
    options.nSamples (1,1) double = 7
    options.nSamplesLargest (1,1) double = 3
    options.seed (1,1) double = 37
    options.errorTolerance (1,1) double = 1e-12
    options.shouldBuild (1,1) logical = true
    options.shouldWriteArtifacts (1,1) logical = true
    options.runId (1,1) string = ""
end

validateOptions(options);

sourceDirectory = fileparts(mfilename('fullpath'));
fftwLibraryPath = fullfile(matlabroot,'bin',computer('arch'),'libmwfftw3.3.dylib');
runId = options.runId;
if strlength(runId) == 0
    timestamp = string(datetime('now',TimeZone='UTC',Format="yyyyMMdd'T'HHmmssSSS'Z'"));
    runId = timestamp + "-" + string(computer('arch')) + "-r" + lower(string(version('-release')));
end

runDirectory = "";
if options.shouldWriteArtifacts
    runDirectory = fullfile(options.outputDirectory,runId);
    if isfolder(runDirectory) || isfile(runDirectory)
        error('FFTWBenchmark:OutputExists','Benchmark output already exists at %s.',runDirectory);
    end
    [didCreate,message] = mkdir(runDirectory);
    if ~didCreate
        error('FFTWBenchmark:OutputCreationFailed','Unable to create %s: %s',runDirectory,message);
    end
end

previousPlanner = string(fftw('planner'));
previousWisdom = fftw('dwisdom');
previousThreads = maxNumCompThreads;
previousRandomState = rng;
stateCleanup = onCleanup(@() restoreMatlabState(previousPlanner,previousWisdom,previousThreads,previousRandomState));

if isnan(options.nThreads)
    maxNumCompThreads('automatic');
    requestedThreads = maxNumCompThreads;
    threadPolicy = "hardwareMaximum";
else
    maxNumCompThreads(options.nThreads);
    requestedThreads = options.nThreads;
    threadPolicy = "explicit";
end
fftw('planner',char(options.planner));
fftw('dwisdom',[]);
observedThreads = maxNumCompThreads;

result = initializeResult(runId,runDirectory,sourceDirectory,fftwLibraryPath,options,threadPolicy,requestedThreads,observedThreads);
activeWorkloadSize = [];

try
    if options.shouldBuild
        clear fftw_dft2
        RealToComplexTransform.makeMexFiles(fftwLibraryPath);
        rehash;
    elseif exist('fftw_dft2','file') ~= 3
        error('FFTWBenchmark:MexMissing','fftw_dft2 is not available and shouldBuild is false.');
    end

    for iSize = 1:size(options.sizes,1)
        activeWorkloadSize = options.sizes(iSize,:);
        nSamples = options.nSamples;
        if iSize == size(options.sizes,1)
            nSamples = options.nSamplesLargest;
        end
        workloadSeed = options.seed + iSize - 1;
        workload = benchmarkWorkload(activeWorkloadSize,workloadSeed,options.nWarmups,nSamples,requestedThreads,options.planner,options.errorTolerance);
        result.workloads(end+1,1) = workload;
        if ~workload.correctness.passed
            error('FFTWBenchmark:CorrectnessFailure','Correctness validation failed for workload [%s].',num2str(activeWorkloadSize));
        end
    end

    result.status = "passed";
    result.completedAtUTC = utcTimestamp();
    if options.shouldWriteArtifacts
        writeResultArtifacts(result,runDirectory);
    end
catch exception
    result.status = "failed";
    result.completedAtUTC = utcTimestamp();
    result.failure = failureRecord(exception,activeWorkloadSize);
    if options.shouldWriteArtifacts
        try
            writeResultArtifacts(result,runDirectory);
        catch artifactException
            warning('FFTWBenchmark:FailureArtifactWriteFailed','Unable to write failure artifacts: %s',artifactException.message);
        end
    end
    rethrow(exception);
end

clear stateCleanup
end

function validateOptions(options)
validateattributes(options.sizes,{'double'},{'2d','ncols',3,'integer','positive'},mfilename,'sizes');
if any(mod(options.sizes(:,2),2) ~= 0)
    error('FFTWBenchmark:OddHalfSpectrumDimension','The current half-y inverse requires even values of Ny.');
end
if ~ismember(options.planner,["estimate","measure","patient","exhaustive"])
    error('FFTWBenchmark:UnknownPlanner','Planner must be estimate, measure, patient, or exhaustive.');
end
if ~isnan(options.nThreads)
    validateattributes(options.nThreads,{'double'},{'integer','positive'},mfilename,'nThreads');
end
validateattributes(options.nWarmups,{'double'},{'integer','nonnegative'},mfilename,'nWarmups');
validateattributes(options.nSamples,{'double'},{'integer','positive'},mfilename,'nSamples');
validateattributes(options.nSamplesLargest,{'double'},{'integer','positive'},mfilename,'nSamplesLargest');
validateattributes(options.seed,{'double'},{'integer','nonnegative'},mfilename,'seed');
validateattributes(options.errorTolerance,{'double'},{'real','finite','nonnegative'},mfilename,'errorTolerance');
if options.shouldWriteArtifacts && strlength(options.outputDirectory) == 0
    error('FFTWBenchmark:MissingOutputDirectory','outputDirectory must be provided when shouldWriteArtifacts is true.');
end
end

function result = initializeResult(runId,runDirectory,sourceDirectory,fftwLibraryPath,options,threadPolicy,requestedThreads,observedThreads)
result.schemaVersion = "1.0.0";
result.status = "running";
result.runId = runId;
result.generatedAtUTC = utcTimestamp();
result.completedAtUTC = "";
result.environment = collectEnvironment(fftwLibraryPath,observedThreads);
result.configuration.sizes = options.sizes;
result.configuration.transformDimensions = [1 2];
result.configuration.halfSpectrumDimension = 2;
result.configuration.inputDistribution = "randn";
result.configuration.randomGenerator = "twister";
result.configuration.seeds = options.seed + (0:size(options.sizes,1)-1);
result.configuration.planner = options.planner;
result.configuration.threadPolicy = threadPolicy;
result.configuration.requestedThreads = requestedThreads;
result.configuration.observedMatlabThreads = observedThreads;
result.configuration.mexThreads = requestedThreads;
result.configuration.nWarmups = options.nWarmups;
result.configuration.nSamples = options.nSamples;
result.configuration.nSamplesLargest = options.nSamplesLargest;
result.configuration.timingBoundary = "transform-call-only";
result.configuration.normalizationIncludedInTiming = false;
result.configuration.operationOrder = "round-robin with a rotating first operation";
result.configuration.errorMetric = "relative infinity norm";
result.configuration.errorTolerance = options.errorTolerance;
result.sources = collectSourceProvenance(sourceDirectory);
result.workloads = repmat(emptyWorkload(),0,1);
result.failure = [];
result.artifacts.directory = runDirectory;
result.artifacts.json = "benchmark.json";
result.artifacts.markdown = "summary.md";
end

function workload = benchmarkWorkload(sz,workloadSeed,nWarmups,nSamples,nThreads,planner,errorTolerance)
rng(workloadSeed,'twister');
x = randn(sz);
dft = RealToComplexTransform(sz,dims=[1 2],nCores=nThreads,planner=char(planner));

fullSpectrumReference = fft(fft(x,sz(1),1),sz(2),2);
halfSpectrumReference = fullSpectrumReference(:,1:dft.complexSize(2),:);
preservingInputBefore = complex(real(halfSpectrumReference),imag(halfSpectrumReference));

matlabForwardOutput = fullSpectrumReference;
matlabInverseOutput = zeros(sz);
mexForwardAllocatingOutput = complex(zeros(dft.complexSize));
mexForwardPreallocatedOutput = complex(zeros(dft.complexSize));
mexInverseAllocatingOutput = zeros(sz);
mexInversePreallocatedOutput = zeros(sz);
mexInverseDestructiveOutput = zeros(sz);
destructiveSpectrum = complex(zeros(dft.complexSize));

operationIds = ["matlab_forward","matlab_inverse","mex_r2c_allocating","mex_r2c_preallocated","mex_c2r_allocating_preserving","mex_c2r_preallocated_preserving","mex_c2r_preallocated_destructive"];
nOperations = numel(operationIds);
sampleTimes = nan(nOperations,nSamples);
nRounds = nWarmups + nSamples;

for iRound = 1:nRounds
    firstOperation = mod(iRound-1,nOperations) + 1;
    operationOrder = [firstOperation:nOperations 1:firstOperation-1];
    for iOperation = operationOrder
        switch iOperation
            case 1
                timer = tic;
                matlabForwardOutput = fft(fft(x,sz(1),1),sz(2),2);
                elapsed = toc(timer);
            case 2
                timer = tic;
                matlabInverseOutput = ifft(ifft(fullSpectrumReference,sz(1),1),sz(2),2,'symmetric');
                elapsed = toc(timer);
            case 3
                timer = tic;
                mexForwardAllocatingOutput = dft.transformForward(x);
                elapsed = toc(timer);
            case 4
                timer = tic;
                mexForwardPreallocatedOutput = dft.transformForwardIntoArray(x,mexForwardPreallocatedOutput);
                elapsed = toc(timer);
            case 5
                timer = tic;
                mexInverseAllocatingOutput = dft.transformBack(halfSpectrumReference);
                elapsed = toc(timer);
            case 6
                timer = tic;
                mexInversePreallocatedOutput = dft.transformBackIntoArray(halfSpectrumReference,mexInversePreallocatedOutput);
                elapsed = toc(timer);
            case 7
                destructiveSpectrum = complex(real(halfSpectrumReference),imag(halfSpectrumReference));
                timer = tic;
                [destructiveSpectrum,mexInverseDestructiveOutput] = dft.transformBackIntoArrayDestructive(destructiveSpectrum,mexInverseDestructiveOutput);
                elapsed = toc(timer);
        end
        if iRound > nWarmups
            sampleTimes(iOperation,iRound-nWarmups) = elapsed;
        end
    end
end

preservingInputUnchanged = isequaln(halfSpectrumReference,preservingInputBefore);
destructiveInputChanged = ~isequaln(destructiveSpectrum,halfSpectrumReference);

storage = storageRecord(x,fullSpectrumReference,halfSpectrumReference,mexInversePreallocatedOutput,[]);
operations = repmat(emptyOperation(),nOperations,1);
operations(1) = operationRecord(operationIds(1),"MATLAB forward","forward","matlab","allocating","not-applicable",sampleTimes(1,:),matlabForwardOutput,fullSpectrumReference,storage,0,true,errorTolerance);
operations(2) = operationRecord(operationIds(2),"MATLAB inverse","inverse","matlab","allocating","not-applicable",sampleTimes(2,:),matlabInverseOutput,x,storage,0,true,errorTolerance);
operations(3) = operationRecord(operationIds(3),"MEX r2c allocating","forward","mex","allocating","not-applicable",sampleTimes(3,:),mexForwardAllocatingOutput,halfSpectrumReference,storage,0,true,errorTolerance);
operations(4) = operationRecord(operationIds(4),"MEX r2c preallocated","forward","mex","preallocated","not-applicable",sampleTimes(4,:),mexForwardPreallocatedOutput,halfSpectrumReference,storage,0,true,errorTolerance);
operations(5) = operationRecord(operationIds(5),"MEX c2r allocating preserving","inverse","mex","allocating","preserving",sampleTimes(5,:),dft.scaleFactor*mexInverseAllocatingOutput,x,storage,0,preservingInputUnchanged,errorTolerance);
operations(6) = operationRecord(operationIds(6),"MEX c2r preallocated preserving","inverse","mex","preallocated","preserving",sampleTimes(6,:),dft.scaleFactor*mexInversePreallocatedOutput,x,storage,0,preservingInputUnchanged,errorTolerance);
operations(7) = operationRecord(operationIds(7),"MEX c2r preallocated destructive","inverse","mex","preallocated","destructive",sampleTimes(7,:),dft.scaleFactor*mexInverseDestructiveOutput,x,storage,storage.mexHalfSpectrumBytes,destructiveInputChanged,errorTolerance);

operations = addAllocationModel(operations,storage);
matlabForwardMedian = operations(1).medianSeconds;
matlabInverseMedian = operations(2).medianSeconds;
for iOperation = 1:nOperations
    if operations(iOperation).direction == "forward"
        operations(iOperation).speedRatioToMatlab = matlabForwardMedian/operations(iOperation).medianSeconds;
    else
        operations(iOperation).speedRatioToMatlab = matlabInverseMedian/operations(iOperation).medianSeconds;
    end
end

workload = emptyWorkload();
workload.size = sz;
workload.seed = workloadSeed;
workload.nWarmups = nWarmups;
workload.nSamples = nSamples;
workload.complexSize = dft.complexSize;
workload.scaleFactor = dft.scaleFactor;
workload.storage = storage;
workload.operations = operations;
workload.correctness.maximumAbsoluteError = max([operations.maximumAbsoluteError]);
workload.correctness.maximumRelativeError = max([operations.maximumRelativeError]);
workload.correctness.preservingInputUnchanged = preservingInputUnchanged;
workload.correctness.destructiveInputChanged = destructiveInputChanged;
workload.correctness.passed = all([operations.correctnessPassed]);
end

function operation = operationRecord(id,label,direction,backend,allocationMode,inputContract,sampleTimes,actual,reference,storage,untimedPreparationBytes,inputContractSatisfied,errorTolerance)
[maximumAbsoluteError,maximumRelativeError] = numericalError(actual,reference);
operation = emptyOperation();
operation.id = id;
operation.label = label;
operation.direction = direction;
operation.backend = backend;
operation.allocationMode = allocationMode;
operation.inputContract = inputContract;
operation.inputContractSatisfied = inputContractSatisfied;
operation.sampleTimesSeconds = sampleTimes;
operation.medianSeconds = median(sampleTimes);
operation.speedRatioToMatlab = NaN;
operation.maximumAbsoluteError = maximumAbsoluteError;
operation.maximumRelativeError = maximumRelativeError;
operation.correctnessPassed = maximumRelativeError <= errorTolerance && inputContractSatisfied;
operation.untimedPreparationBytesPerCall = untimedPreparationBytes;
if backend == "mex"
    operation.persistentTransformBytes = storage.totalPersistentTransformBytes;
end
end

function operations = addAllocationModel(operations,storage)
operations(1).minimumKnownAllocationCount = 2;
operations(1).minimumKnownAllocatedBytesPerCall = 2*storage.matlabFullSpectrumBytes;
operations(1).unresolvedAllocationSources = "MATLAB internal FFT work buffers";

operations(2).minimumKnownAllocationCount = 2;
operations(2).minimumKnownAllocatedBytesPerCall = storage.matlabFullSpectrumBytes + storage.realOutputBytes;
operations(2).unresolvedAllocationSources = "MATLAB internal FFT work buffers";

operations(3).minimumKnownAllocationCount = 1;
operations(3).minimumKnownAllocatedBytesPerCall = storage.mexHalfSpectrumBytes;

operations(4).callerPreallocatedBytes = storage.mexHalfSpectrumBytes;
operations(4).unresolvedAllocationSources = "MATLAB complex(...) wrapper and copy-on-write behavior";

operations(5).minimumKnownAllocationCount = 1;
operations(5).minimumKnownAllocatedBytesPerCall = storage.realOutputBytes;
operations(5).minimumKnownCopiedBytesPerCall = storage.mexHalfSpectrumBytes;
operations(5).unresolvedAllocationSources = "MATLAB complex(...) wrapper and copy-on-write behavior";

operations(6).callerPreallocatedBytes = storage.realOutputBytes;
operations(6).minimumKnownCopiedBytesPerCall = storage.mexHalfSpectrumBytes;
operations(6).unresolvedAllocationSources = "MATLAB complex(...) wrapper and copy-on-write behavior";

operations(7).callerPreallocatedBytes = storage.mexHalfSpectrumBytes + storage.realOutputBytes;
operations(7).unresolvedAllocationSources = "MATLAB complex(...) wrapper and copy-on-write behavior";
end

function storage = storageRecord(x,fullSpectrum,halfSpectrum,realOutput,wrapperScratch)
storage.realInputBytes = numel(x)*8;
storage.matlabFullSpectrumBytes = numel(fullSpectrum)*16;
storage.mexHalfSpectrumBytes = numel(halfSpectrum)*16;
storage.realOutputBytes = numel(realOutput)*8;
storage.wrapperScratchBytes = numel(wrapperScratch)*16;
storage.mexInternalScratchBytes = storage.mexHalfSpectrumBytes;
storage.totalPersistentTransformBytes = storage.wrapperScratchBytes + storage.mexInternalScratchBytes;
storage.halfToFullSpectrumStorageRatio = storage.mexHalfSpectrumBytes/storage.matlabFullSpectrumBytes;
end

function [maximumAbsoluteError,maximumRelativeError] = numericalError(actual,reference)
maximumAbsoluteError = max(abs(actual-reference),[],'all');
referenceMagnitude = max(abs(reference),[],'all');
if referenceMagnitude == 0
    referenceMagnitude = 1;
end
maximumRelativeError = maximumAbsoluteError/referenceMagnitude;
end

function environment = collectEnvironment(fftwLibraryPath,observedThreads)
environment.matlabVersion = string(version);
environment.matlabRelease = string(version('-release'));
environment.architecture = string(computer('arch'));
environment.mexExtension = string(mexext);
environment.operatingSystem = string(system_dependent('getos'));
environment.processor = "unknown";
environment.machineModel = "unknown";
environment.machineName = "unknown";
environment.physicalMemory = "unknown";

[hardwareStatus,hardwareText] = system('system_profiler SPHardwareDataType -json');
if hardwareStatus == 0
    hardwareData = jsondecode(hardwareText);
    hardware = hardwareData.SPHardwareDataType(1);
    environment.processor = optionalStructString(hardware,'chip_type',environment.processor);
    environment.machineModel = optionalStructString(hardware,'machine_model',environment.machineModel);
    environment.machineName = optionalStructString(hardware,'machine_name',environment.machineName);
    environment.physicalMemory = optionalStructString(hardware,'physical_memory',environment.physicalMemory);
end

environment.observedComputationalThreads = observedThreads;
environment.fftwLibrary = string(fftwLibraryPath);
wisdom = fftw('dwisdom');
fftwVersion = regexp(wisdom,'fftw-[0-9.]+','match','once');
if isempty(fftwVersion)
    environment.fftwVersion = "unknown";
else
    environment.fftwVersion = string(fftwVersion);
end

environment.mexCompiler = "unknown";
environment.mexCompilerVersion = "unknown";
try
    compiler = mex.getCompilerConfigurations('C++','Selected');
    if ~isempty(compiler)
        environment.mexCompiler = string(compiler(1).Name);
        environment.mexCompilerVersion = string(compiler(1).Version);
    end
catch
end
end

function value = optionalStructString(inputStruct,fieldName,defaultValue)
if isfield(inputStruct,fieldName)
    value = string(inputStruct.(fieldName));
else
    value = defaultValue;
end
end

function sources = collectSourceProvenance(sourceDirectory)
sourcePaths = [fullfile(sourceDirectory,"runFFTWFeasibilityBaseline.m"), fullfile(sourceDirectory,"RealToComplexTransform.m"), fullfile(sourceDirectory,"fftw_dft2.cpp")];
sourceNames = ["runFFTWFeasibilityBaseline.m","RealToComplexTransform.m","fftw_dft2.cpp"];
files = repmat(struct('path',"",'sha256',""),numel(sourcePaths),1);
for iFile = 1:numel(sourcePaths)
    files(iFile).path = sourceNames(iFile);
    files(iFile).sha256 = fileSHA256(sourcePaths(iFile));
end
sources.files = files;

repositoryRoot = fileparts(fileparts(fileparts(sourceDirectory)));
[commitStatus,commitText] = system(sprintf('git -C "%s" rev-parse HEAD',repositoryRoot));
if commitStatus == 0
    sources.repositoryCommit = string(strtrim(commitText));
else
    sources.repositoryCommit = "unknown";
end
[dirtyStatus,dirtyText] = system(sprintf('git -C "%s" status --porcelain --untracked-files=all',repositoryRoot));
sources.repositoryDirty = dirtyStatus ~= 0 || strlength(strtrim(string(dirtyText))) > 0;
end

function hash = fileSHA256(path)
digest = java.security.MessageDigest.getInstance('SHA-256');
bytes = unicode2native(fileread(path),'UTF-8');
digest.update(bytes);
hashBytes = typecast(digest.digest(),'uint8');
hash = string(lower(reshape(dec2hex(hashBytes,2).',1,[])));
end

function writeResultArtifacts(result,runDirectory)
serializableResult = result;
serializableResult.workloads = num2cell(result.workloads);
jsonText = jsonencode(serializableResult,PrettyPrint=true);
writeTextFile(fullfile(runDirectory,"benchmark.json"),string(jsonText));
writeTextFile(fullfile(runDirectory,"summary.md"),createMarkdownSummary(result));
end

function summary = createMarkdownSummary(result)
lines = strings(0,1);
lines(end+1) = "# FFTW r2c feasibility baseline";
lines(end+1) = "";
lines(end+1) = "Status: **" + upper(result.status) + "**";
lines(end+1) = "";
lines(end+1) = "This report is a descriptive Stage 1 baseline. It does not make a GO or NO-GO decision.";
lines(end+1) = "";
lines(end+1) = "## Environment";
lines(end+1) = "";
lines(end+1) = "| Field | Value |";
lines(end+1) = "|---|---|";
lines(end+1) = markdownRow("Run",result.runId);
lines(end+1) = markdownRow("Generated (UTC)",result.generatedAtUTC);
lines(end+1) = markdownRow("MATLAB",result.environment.matlabVersion);
lines(end+1) = markdownRow("Architecture",result.environment.architecture);
lines(end+1) = markdownRow("Operating system",result.environment.operatingSystem);
lines(end+1) = markdownRow("Processor",result.environment.processor);
lines(end+1) = markdownRow("Machine",result.environment.machineName + " " + result.environment.machineModel);
lines(end+1) = markdownRow("Physical memory",result.environment.physicalMemory);
lines(end+1) = markdownRow("FFTW",result.environment.fftwVersion + " via " + result.environment.fftwLibrary);
lines(end+1) = markdownRow("Planner",result.configuration.planner);
lines(end+1) = markdownRow("Threads",string(result.configuration.observedMatlabThreads));
lines(end+1) = markdownRow("Warmups",string(result.configuration.nWarmups));
lines(end+1) = markdownRow("Samples",sprintf('%d (%d for largest)',result.configuration.nSamples,result.configuration.nSamplesLargest));

if ~isempty(result.workloads)
    lines(end+1) = "";
    lines(end+1) = "## Median transform-call time";
    lines(end+1) = "";
    lines(end+1) = "MEX cells show seconds followed by speed relative to the corresponding MATLAB transform. Ratios above 1 are faster than MATLAB.";
    lines(end+1) = "";
    lines(end+1) = "| Size | MATLAB forward | MEX r2c alloc | MEX r2c prealloc | MATLAB inverse | MEX c2r alloc | MEX c2r prealloc | MEX c2r destructive |";
    lines(end+1) = "|---|---:|---:|---:|---:|---:|---:|---:|";
    for iWorkload = 1:numel(result.workloads)
        operations = result.workloads(iWorkload).operations;
        lines(end+1) = "| " + formatSize(result.workloads(iWorkload).size) + " | " + formatTiming(operations(1),false) + " | " + formatTiming(operations(3),true) + " | " + formatTiming(operations(4),true) + " | " + formatTiming(operations(2),false) + " | " + formatTiming(operations(5),true) + " | " + formatTiming(operations(6),true) + " | " + formatTiming(operations(7),true) + " |"; %#ok<AGROW>
    end

    lines(end+1) = "";
    lines(end+1) = "## Correctness";
    lines(end+1) = "";
    lines(end+1) = "The pass threshold is a relative infinity error of `" + string(result.configuration.errorTolerance) + "`.";
    lines(end+1) = "";
    lines(end+1) = "| Size | Maximum absolute error | Maximum relative error | Preserving input unchanged | Destructive input changed | Pass |";
    lines(end+1) = "|---|---:|---:|:---:|:---:|:---:|";
    for iWorkload = 1:numel(result.workloads)
        correctness = result.workloads(iWorkload).correctness;
        lines(end+1) = sprintf('| %s | %.6g | %.6g | %s | %s | %s |',formatSize(result.workloads(iWorkload).size),correctness.maximumAbsoluteError,correctness.maximumRelativeError,yesNo(correctness.preservingInputUnchanged),yesNo(correctness.destructiveInputChanged),yesNo(correctness.passed)); %#ok<AGROW>
    end

    lines(end+1) = "";
    lines(end+1) = "## Array storage";
    lines(end+1) = "";
    lines(end+1) = "| Size | Real array (MiB) | MATLAB full spectrum (MiB) | MEX half spectrum (MiB) | Half/full ratio | Persistent transform buffers (MiB) |";
    lines(end+1) = "|---|---:|---:|---:|---:|---:|";
    for iWorkload = 1:numel(result.workloads)
        storage = result.workloads(iWorkload).storage;
        lines(end+1) = sprintf('| %s | %.3f | %.3f | %.3f | %.6f | %.3f |',formatSize(result.workloads(iWorkload).size),toMiB(storage.realInputBytes),toMiB(storage.matlabFullSpectrumBytes),toMiB(storage.mexHalfSpectrumBytes),storage.halfToFullSpectrumStorageRatio,toMiB(storage.totalPersistentTransformBytes)); %#ok<AGROW>
    end
end

lines(end+1) = "";
lines(end+1) = "## Allocation model";
lines(end+1) = "";
lines(end+1) = "Counts are source-grounded minimums. MATLAB internal FFT work buffers and possible `complex(...)` or copy-on-write allocations remain unresolved for issue #39.";
lines(end+1) = "";
lines(end+1) = "| Operation | Minimum known timed allocation | Minimum known timed copy | Caller preallocation |";
lines(end+1) = "|---|---|---|---|";
lines(end+1) = "| MATLAB forward | Two full-spectrum arrays | Unknown internal work | None |";
lines(end+1) = "| MATLAB inverse | One full-spectrum and one real array | Unknown internal work | None |";
lines(end+1) = "| MEX r2c allocating | One half-spectrum output | None known | None |";
lines(end+1) = "| MEX r2c preallocated | None known | Copy-on-write unresolved | Half-spectrum output |";
lines(end+1) = "| MEX c2r allocating preserving | One real output | One half-spectrum `memcpy` | None |";
lines(end+1) = "| MEX c2r preallocated preserving | None known | One half-spectrum `memcpy` | Real output |";
lines(end+1) = "| MEX c2r preallocated destructive | None known | Copy-on-write unresolved | Half-spectrum input and real output |";

if result.status == "failed" && ~isempty(result.failure)
    lines(end+1) = "";
    lines(end+1) = "## Failure";
    lines(end+1) = "";
    lines(end+1) = "- Identifier: `" + result.failure.identifier + "`";
    lines(end+1) = "- Message: " + escapeMarkdown(result.failure.message);
    if ~isempty(result.failure.activeWorkloadSize)
        lines(end+1) = "- Active workload: `" + formatSize(result.failure.activeWorkloadSize) + "`";
    end
end

summary = strjoin(lines,newline);
end

function row = markdownRow(field,value)
row = "| " + escapeMarkdown(string(field)) + " | " + escapeMarkdown(string(value)) + " |";
end

function value = escapeMarkdown(value)
value = replace(value,"|","\|");
value = replace(value,newline,"<br>");
end

function value = formatTiming(operation,includeRatio)
if includeRatio
    value = sprintf('%.6g (%.3fx)',operation.medianSeconds,operation.speedRatioToMatlab);
else
    value = sprintf('%.6g',operation.medianSeconds);
end
end

function value = formatSize(sz)
value = strtrim(sprintf('%d x %d x %d',sz));
end

function value = yesNo(tf)
if tf
    value = 'yes';
else
    value = 'no';
end
end

function value = toMiB(bytes)
value = bytes/(1024^2);
end

function failure = failureRecord(exception,activeWorkloadSize)
failure.identifier = string(exception.identifier);
failure.message = string(exception.message);
failure.activeWorkloadSize = activeWorkloadSize;
stack = exception.stack;
failure.stack = repmat(struct('name',"",'line',0),numel(stack),1);
for iFrame = 1:numel(stack)
    failure.stack(iFrame).name = string(stack(iFrame).name);
    failure.stack(iFrame).line = stack(iFrame).line;
end
end

function writeTextFile(path,text)
[fileId,message] = fopen(path,'w');
if fileId < 0
    error('FFTWBenchmark:ArtifactOpenFailed','Unable to open %s: %s',path,message);
end
fileCleanup = onCleanup(@() fclose(fileId));
fprintf(fileId,'%s\n',text);
clear fileCleanup
end

function restoreMatlabState(previousPlanner,previousWisdom,previousThreads,previousRandomState)
try
    fftw('dwisdom',[]);
    fftw('dwisdom',previousWisdom);
catch
end
try
    fftw('planner',char(previousPlanner));
catch
end
try
    maxNumCompThreads(previousThreads);
catch
end
try
    rng(previousRandomState);
catch
end
end

function timestamp = utcTimestamp()
timestamp = string(datetime('now',TimeZone='UTC',Format="yyyy-MM-dd'T'HH:mm:ss.SSS'Z'"));
end

function operation = emptyOperation()
operation.id = "";
operation.label = "";
operation.direction = "";
operation.backend = "";
operation.allocationMode = "";
operation.inputContract = "";
operation.inputContractSatisfied = true;
operation.sampleTimesSeconds = [];
operation.medianSeconds = NaN;
operation.speedRatioToMatlab = NaN;
operation.maximumAbsoluteError = NaN;
operation.maximumRelativeError = NaN;
operation.correctnessPassed = false;
operation.minimumKnownAllocationCount = 0;
operation.minimumKnownAllocatedBytesPerCall = 0;
operation.minimumKnownCopiedBytesPerCall = 0;
operation.callerPreallocatedBytes = 0;
operation.untimedPreparationBytesPerCall = 0;
operation.persistentTransformBytes = 0;
operation.unresolvedAllocationSources = strings(0,1);
end

function workload = emptyWorkload()
workload.size = [];
workload.seed = [];
workload.nWarmups = [];
workload.nSamples = [];
workload.complexSize = [];
workload.scaleFactor = [];
workload.storage = struct;
workload.operations = repmat(emptyOperation(),0,1);
workload.correctness.maximumAbsoluteError = [];
workload.correctness.maximumRelativeError = [];
workload.correctness.preservingInputUnchanged = [];
workload.correctness.destructiveInputChanged = [];
workload.correctness.passed = false;
end
