function result = runFFTWEngineLayoutBenchmark(options)
% Compare FFT engines, half-spectrum layouts, planning, and threading.
%
% This experimental issue #38 entry point screens the complete requested
% configuration matrix and then remeasures each workload winner alongside
% MATLAB. It reports feasibility thresholds without making a GO or NO-GO
% decision for the milestone.
arguments
    options.outputDirectory (1,1) string = fullfile(fileparts(mfilename('fullpath')),"results","issue38")
    options.sizes (:,3) double = [128 128 128; 256 128 40; 256 256 64; 1024 1024 30]
    options.gateSizes (:,3) double = [256 256 64; 1024 1024 30]
    options.transformOrders (:,2) double = [1 2; 2 1]
    options.fftwStrategies (1,:) string = ["guru-rank2","staged-r2c-c2c"]
    options.vdspStrategies (1,:) string = ["vdsp-2d","vdsp-staged"]
    options.planners (1,:) string = ["estimate","measure","patient","exhaustive"]
    options.threadCounts (1,:) double = [1 2 4 6 8 12 18]
    options.alignmentModes (1,:) string = ["matched","unaligned"]
    options.plannerTimeLimitSeconds (1,1) double = 10
    options.nScreeningWarmups (1,1) double = 1
    options.nScreeningSamples (1,1) double = 3
    options.nWarmups (1,1) double = 2
    options.nSamples (1,1) double = 7
    options.nSamplesLargest (1,1) double = 3
    options.seed (1,1) double = 38
    options.errorTolerance (1,1) double = 1e-12
    options.rawForwardSpeedThreshold (1,1) double = 1.25
    options.totalForwardSpeedThreshold (1,1) double = 1.10
    options.inverseSpeedThreshold (1,1) double = 0.95
    options.fallbackPolicy (1,1) string = "threshold"
    options.shouldBuild (1,1) logical = true
    options.shouldBuildNative (1,1) logical = true
    options.shouldWriteArtifacts (1,1) logical = true
    options.runId (1,1) string = ""
end

validateOptions(options);
sourceDirectory = fileparts(mfilename('fullpath'));
runId = options.runId;
if strlength(runId) == 0
    timestamp = string(datetime('now',TimeZone='UTC',Format="yyyyMMdd'T'HHmmssSSS'Z'"));
    runId = timestamp + "-" + string(computer('arch')) + "-r" + lower(string(version('-release')));
end
runDirectory = prepareRunDirectory(options.outputDirectory,runId,options.shouldWriteArtifacts);

previousPlanner = string(fftw('planner'));
previousWisdom = fftw('dwisdom');
previousThreads = maxNumCompThreads;
previousRandomState = rng;
stateCleanup = onCleanup(@() restoreMatlabState(previousPlanner,previousWisdom,previousThreads,previousRandomState));
maxNumCompThreads('automatic');
hardwareThreads = maxNumCompThreads;
requestedThreads = unique(min(options.threadCounts,hardwareThreads),'stable');
requestedThreads = requestedThreads(requestedThreads >= 1);
fftw('planner','measure');
fftw('dwisdom',[]);

result = initializeResult(runId,runDirectory,sourceDirectory,options,requestedThreads,hardwareThreads);
activeStage = "build";
try
    if options.shouldBuild
        result.build = buildFFTWEngineBenchmarkMex;
    else
        requireModule("fftw_engine_benchmark_bundled");
        result.build = existingBuildRecord();
    end
    result.engines = collectInitialEngineRecords(result.build);
    [bundledVersion,bundledLibrary] = fftw_engine_benchmark_bundled('info');
    result.engines(1).version = string(bundledVersion);
    result.engines(1).library = string(bundledLibrary);
    result.alignmentValidation(end+1,1) = alignmentValidationRecord("bundled-fftw","fftw_engine_benchmark_bundled");

    activeStage = "bundled FFTW screening";
    bundledCandidates = screenFFTWEngine("bundled-fftw","fftw_engine_benchmark_bundled",options,requestedThreads);
    result.candidates = [result.candidates; bundledCandidates];
    bundledFinalists = measureEngineFinalists(result.candidates,"bundled-fftw",options);
    result.bundledAssessment = assessFinalists(bundledFinalists,options);

    shouldRunFallbacks = options.fallbackPolicy == "always" || (options.fallbackPolicy == "threshold" && ~result.bundledAssessment.forwardGatePassed);
    result.configuration.fallbackTriggered = shouldRunFallbacks;
    if shouldRunFallbacks
        activeStage = "native FFTW build and screening";
        canRunNative = false;
        if options.shouldBuildNative
            try
                nativeBuild = buildFFTWEngineBenchmarkMex(shouldBuildNative=true,shouldBuildBundledAndVDSP=false);
                result.build.native = nativeBuild.native;
                canRunNative = true;
            catch exception
                result.engines(end+1,1) = unavailableEngine("native-fftw",exception);
            end
        end
        if canRunNative
            [nativeVersion,nativeLibrary] = fftw_engine_benchmark_native('info');
            if startsWith(string(nativeLibrary),string(matlabroot))
                result.engines(end+1,1) = unavailableEngineMessage("native-fftw","The native module resolved FFTW symbols to MATLAB's bundled library and was rejected.");
            else
                result.engines(end+1,1) = availableEngine("native-fftw","fftw_engine_benchmark_native",string(nativeVersion),string(nativeLibrary));
                result.alignmentValidation(end+1,1) = alignmentValidationRecord("native-fftw","fftw_engine_benchmark_native");
                nativeCandidates = screenFFTWEngine("native-fftw","fftw_engine_benchmark_native",options,requestedThreads);
                result.candidates = [result.candidates; nativeCandidates];
            end
        elseif ~any([result.engines.id] == "native-fftw")
            result.engines(end+1,1) = unavailableEngineMessage("native-fftw","The native FFTW module was not built.");
        end

        activeStage = "Accelerate/vDSP screening";
        if ismac && exist('vdsp_engine_benchmark','file') == 3
            [vdspVersion,vdspLibrary] = vdsp_engine_benchmark('info');
            result.engines(end+1,1) = availableEngine("accelerate-vdsp","vdsp_engine_benchmark",string(vdspVersion),string(vdspLibrary));
            vdspCandidates = screenVDSPEngine(options);
            result.candidates = [result.candidates; vdspCandidates];
        else
            result.engines(end+1,1) = unavailableEngineMessage("accelerate-vdsp","Accelerate/vDSP is unavailable on this platform.");
        end
    end

    activeStage = "winner measurement";
    result.workloads = measureGlobalFinalists(result.candidates,options);
    result.thresholdSummary = assessWorkloads(result.workloads,options);
    result.status = "passed";
    result.completedAtUTC = utcTimestamp();
    if options.shouldWriteArtifacts
        writeResultArtifacts(result,runDirectory);
    end
catch exception
    result.status = "failed";
    result.completedAtUTC = utcTimestamp();
    result.failure = failureRecord(exception,activeStage);
    if options.shouldWriteArtifacts
        try
            writeResultArtifacts(result,runDirectory);
        catch artifactException
            warning('FFTWEngineBenchmark:FailureArtifactWriteFailed','Unable to write issue #38 failure artifacts: %s',artifactException.message);
        end
    end
    rethrow(exception);
end
clear stateCleanup
end

function validateOptions(options)
validateattributes(options.sizes,{'double'},{'2d','ncols',3,'integer','positive'},mfilename,'sizes');
validateattributes(options.gateSizes,{'double'},{'2d','ncols',3,'integer','positive'},mfilename,'gateSizes');
validateattributes(options.transformOrders,{'double'},{'2d','ncols',2,'integer','positive'},mfilename,'transformOrders');
if any(options.transformOrders(:,1) == options.transformOrders(:,2)) || any(options.transformOrders > 3,'all')
    error('FFTWEngineBenchmark:InvalidTransformOrders','Each transform order must contain two distinct dimensions between 1 and 3.');
end
if ~all(ismember(options.fftwStrategies,["guru-rank2","staged-r2c-c2c"]))
    error('FFTWEngineBenchmark:UnknownFFTWStrategy','Unknown FFTW strategy requested.');
end
if ~all(ismember(options.vdspStrategies,["vdsp-2d","vdsp-staged"]))
    error('FFTWEngineBenchmark:UnknownVDSPStrategy','Unknown vDSP strategy requested.');
end
if ~all(ismember(options.planners,["estimate","measure","patient","exhaustive"]))
    error('FFTWEngineBenchmark:UnknownPlanner','Unknown FFTW planner requested.');
end
if ~all(ismember(options.alignmentModes,["matched","unaligned"]))
    error('FFTWEngineBenchmark:UnknownAlignmentMode','Alignment modes must be matched or unaligned.');
end
if ~ismember(options.fallbackPolicy,["threshold","always","never"])
    error('FFTWEngineBenchmark:UnknownFallbackPolicy','fallbackPolicy must be threshold, always, or never.');
end
validateattributes(options.threadCounts,{'double'},{'vector','integer','positive'},mfilename,'threadCounts');
validateattributes(options.plannerTimeLimitSeconds,{'double'},{'scalar','real','finite','positive'},mfilename,'plannerTimeLimitSeconds');
validateattributes(options.nScreeningWarmups,{'double'},{'scalar','integer','nonnegative'},mfilename,'nScreeningWarmups');
validateattributes(options.nScreeningSamples,{'double'},{'scalar','integer','positive'},mfilename,'nScreeningSamples');
validateattributes(options.nWarmups,{'double'},{'scalar','integer','nonnegative'},mfilename,'nWarmups');
validateattributes(options.nSamples,{'double'},{'scalar','integer','positive'},mfilename,'nSamples');
validateattributes(options.nSamplesLargest,{'double'},{'scalar','integer','positive'},mfilename,'nSamplesLargest');
validateattributes(options.errorTolerance,{'double'},{'scalar','real','finite','nonnegative'},mfilename,'errorTolerance');
for iOrder = 1:size(options.transformOrders,1)
    compressedDimension = options.transformOrders(iOrder,2);
    if any(mod(options.sizes(:,compressedDimension),2) ~= 0)
        error('FFTWEngineBenchmark:OddCompressedDimension','Every compressed transform dimension must have even length.');
    end
end
end

function runDirectory = prepareRunDirectory(outputDirectory,runId,shouldWrite)
runDirectory = "";
if ~shouldWrite
    return
end
runDirectory = fullfile(outputDirectory,runId);
if isfolder(runDirectory) || isfile(runDirectory)
    error('FFTWEngineBenchmark:OutputExists','Benchmark output already exists at %s.',runDirectory);
end
[didCreate,message] = mkdir(runDirectory);
if ~didCreate
    error('FFTWEngineBenchmark:OutputCreationFailed','Unable to create %s: %s',runDirectory,message);
end
end

function result = initializeResult(runId,runDirectory,sourceDirectory,options,requestedThreads,hardwareThreads)
result.schemaVersion = "2.0.0";
result.status = "running";
result.runId = runId;
result.generatedAtUTC = utcTimestamp();
result.completedAtUTC = "";
result.environment = collectEnvironment(hardwareThreads);
result.configuration.sizes = options.sizes;
result.configuration.gateSizes = options.gateSizes;
result.configuration.transformOrders = options.transformOrders;
result.configuration.layoutNames = arrayfun(@layoutName,options.transformOrders(:,2));
result.configuration.fftwStrategies = options.fftwStrategies;
result.configuration.vdspStrategies = options.vdspStrategies;
result.configuration.planners = options.planners;
result.configuration.requestedThreadCounts = options.threadCounts;
result.configuration.observedThreadCounts = requestedThreads;
result.configuration.alignmentModes = options.alignmentModes;
result.configuration.plannerTimeLimitSeconds = options.plannerTimeLimitSeconds;
result.configuration.nScreeningWarmups = options.nScreeningWarmups;
result.configuration.nScreeningSamples = options.nScreeningSamples;
result.configuration.nWarmups = options.nWarmups;
result.configuration.nSamples = options.nSamples;
result.configuration.nSamplesLargest = options.nSamplesLargest;
result.configuration.errorTolerance = options.errorTolerance;
result.configuration.rawForwardSpeedThreshold = options.rawForwardSpeedThreshold;
result.configuration.totalForwardSpeedThreshold = options.totalForwardSpeedThreshold;
result.configuration.inverseSpeedThreshold = options.inverseSpeedThreshold;
result.configuration.fallbackPolicy = options.fallbackPolicy;
result.configuration.fallbackTriggered = false;
result.configuration.ranking = "lowest complete forward MEX median; raw pipeline median breaks ties";
result.configuration.rawTimingBoundary = "engine-required pipeline excluding MATLAB allocation and dispatch";
result.configuration.finalTimingOrder = "round-robin with rotating first operation";
result.sources = collectSourceProvenance(sourceDirectory);
result.build = struct;
result.engines = repmat(emptyEngine(),0,1);
result.alignmentValidation = repmat(emptyAlignmentValidation(),0,1);
result.candidates = repmat(emptyCandidate(),0,1);
result.bundledAssessment = emptyAssessment();
result.workloads = repmat(emptyWorkload(),0,1);
result.thresholdSummary = emptyAssessment();
result.failure = [];
result.artifacts.directory = runDirectory;
result.artifacts.json = "engine-layout-benchmark.json";
result.artifacts.markdown = "summary.md";
end

function build = existingBuildRecord()
build.bundled.module = "fftw_engine_benchmark_bundled";
build.bundled.status = "existing";
build.bundled.library = fullfile(matlabroot,"bin",computer('arch'),"libmwfftw3.3.dylib");
build.vdsp.module = "vdsp_engine_benchmark";
build.vdsp.status = string(exist('vdsp_engine_benchmark','file') == 3);
build.native.module = "fftw_engine_benchmark_native";
build.native.status = string(exist('fftw_engine_benchmark_native','file') == 3);
end

function requireModule(module)
if exist(module,'file') ~= 3
    error('FFTWEngineBenchmark:MexMissing','%s is unavailable and shouldBuild is false.',module);
end
end

function engines = collectInitialEngineRecords(build)
engines = availableEngine("bundled-fftw","fftw_engine_benchmark_bundled","unknown",string(build.bundled.library));
end

function candidates = screenFFTWEngine(engine,module,options,threadCounts)
candidates = repmat(emptyCandidate(),0,1);
for iSize = 1:size(options.sizes,1)
    sz = options.sizes(iSize,:);
    rng(options.seed+iSize-1,'twister');
    x = randn(sz);
    for iOrder = 1:size(options.transformOrders,1)
        transformOrder = options.transformOrders(iOrder,:);
        [referenceSpectrum,complexSize] = spectrumReference(x,transformOrder);
        for strategy = options.fftwStrategies
            for planner = options.planners
                for nThreads = threadCounts
                    for alignmentMode = options.alignmentModes
                        feval(module,'forgetWisdom');
                        configuration = candidateConfiguration(engine,module,sz,iSize,transformOrder,strategy,planner,nThreads,alignmentMode);
                        candidate = screenCandidate(configuration,x,referenceSpectrum,complexSize,options);
                        candidates(end+1,1) = candidate; %#ok<AGROW>
                    end
                end
            end
        end
    end
end
end

function candidates = screenVDSPEngine(options)
candidates = repmat(emptyCandidate(),0,1);
for iSize = 1:size(options.sizes,1)
    sz = options.sizes(iSize,:);
    rng(options.seed+iSize-1,'twister');
    x = randn(sz);
    for iOrder = 1:size(options.transformOrders,1)
        transformOrder = options.transformOrders(iOrder,:);
        [referenceSpectrum,complexSize] = spectrumReference(x,transformOrder);
        for strategy = options.vdspStrategies
            configuration = candidateConfiguration("accelerate-vdsp","vdsp_engine_benchmark",sz,iSize,transformOrder,strategy,"not-applicable",0,"16-byte-observed");
            candidate = screenCandidate(configuration,x,referenceSpectrum,complexSize,options);
            candidates(end+1,1) = candidate; %#ok<AGROW>
        end
    end
end
end

function configuration = candidateConfiguration(engine,module,sz,workloadIndex,transformOrder,strategy,planner,nThreads,alignmentMode)
configuration.engine = engine;
configuration.module = module;
configuration.size = sz;
configuration.workloadIndex = workloadIndex;
configuration.transformOrder = transformOrder;
configuration.layout = layoutName(transformOrder(2));
configuration.strategy = strategy;
configuration.planner = planner;
configuration.threads = nThreads;
configuration.alignmentMode = alignmentMode;
end

function candidate = screenCandidate(configuration,x,referenceSpectrum,complexSize,options)
candidate = emptyCandidate();
candidate.engine = configuration.engine;
candidate.module = configuration.module;
candidate.size = configuration.size;
candidate.workloadIndex = configuration.workloadIndex;
candidate.transformOrder = configuration.transformOrder;
candidate.layout = configuration.layout;
candidate.strategy = configuration.strategy;
candidate.planner = configuration.planner;
candidate.threads = configuration.threads;
candidate.alignmentMode = configuration.alignmentMode;
candidate.storage = storageRecord(configuration.size,complexSize);
try
    [plan,planRecord] = createPlan(configuration,x,complexSize,options.plannerTimeLimitSeconds);
    candidate.plan = planRecord;
    planCleanup = onCleanup(@() feval(configuration.module,'free',plan));
    candidate.screening = measureConfiguration(configuration,plan,x,referenceSpectrum,complexSize,options.nScreeningWarmups,options.nScreeningSamples,false);
    candidate.maximumRelativeError = candidate.screening.maximumRelativeError;
    candidate.destructiveInputChanged = candidate.screening.destructiveInputChanged;
    candidate.status = "passed";
    if candidate.maximumRelativeError > options.errorTolerance || ~candidate.destructiveInputChanged || ~isfinite(candidate.screening.forward.totalMexMedianSeconds) || ~isfinite(candidate.screening.inverse.totalMexMedianSeconds)
        candidate.status = "invalid";
        candidate.failure.identifier = "FFTWEngineBenchmark:CandidateValidationFailed";
        candidate.failure.message = "Correctness, destructive semantics, or a required timing path failed.";
    end
    clear planCleanup
catch exception
    candidate.status = "failed";
    candidate.failure = exceptionRecord(exception);
end
end

function [plan,record] = createPlan(configuration,x,complexSize,timeLimit)
spectrumTemplate = complex(zeros(complexSize));
if startsWith(configuration.engine,"accelerate")
    [plan,reportedSize,scaleFactor,planningSeconds,threadPolicy] = feval(configuration.module,'create',configuration.size,configuration.transformOrder,char(configuration.strategy));
    record.planningLimitReached = false;
    record.inputAlignmentClass = pointerAlignment(x);
    record.outputAlignmentClass = pointerAlignment(spectrumTemplate);
    record.observedThreads = threadPolicy;
else
    flags = plannerFlags(configuration.planner);
    [plan,reportedSize,scaleFactor,planningSeconds,limitReached,inputAlignment,outputAlignment,observedThreads] = feval(configuration.module,'create',configuration.size,configuration.transformOrder,char(configuration.strategy),configuration.threads,flags,char(configuration.alignmentMode),timeLimit,x,spectrumTemplate);
    record.planningLimitReached = logical(limitReached);
    record.inputAlignmentClass = inputAlignment;
    record.outputAlignmentClass = outputAlignment;
    record.observedThreads = string(observedThreads);
end
if ~isequal(double(reportedSize),double(complexSize))
    feval(configuration.module,'free',plan);
    error('FFTWEngineBenchmark:ReportedSizeMismatch','The engine reported an unexpected half-spectrum size.');
end
record.planningSeconds = planningSeconds;
record.scaleFactor = scaleFactor;
end

function flags = plannerFlags(planner)
switch planner
    case "estimate"
        flags = 64;
    case "measure"
        flags = 0;
    case "patient"
        flags = 32;
    case "exhaustive"
        flags = 8;
    otherwise
        error('FFTWEngineBenchmark:UnknownPlanner','Unknown planner %s.',planner);
end
end

function measurement = measureConfiguration(configuration,plan,x,referenceSpectrum,complexSize,nWarmups,nSamples,includeMatlab)
operationIds = ["engine-forward-allocating","engine-forward-preallocated","engine-inverse-destructive"];
if includeMatlab
    operationIds = ["matlab-forward","matlab-inverse",operationIds];
end
nOperations = numel(operationIds);
totalSamples = nan(nOperations,nSamples);
kernelSamples = nan(nOperations,nSamples);
pipelineSamples = nan(nOperations,nSamples);
inputAlignmentSamples = nan(nOperations,nSamples);
outputAlignmentSamples = nan(nOperations,nSamples);
operationAvailable = true(1,nOperations);
operationFailures = strings(1,nOperations);
preallocatedSpectrum = complex(zeros(complexSize));
realOutput = zeros(configuration.size);
allocatedOutput = complex(zeros(complexSize));
inverseOutput = zeros(configuration.size);
destroyedSpectrum = complex(referenceSpectrum);
matlabForwardOutput = fftAlong(x,configuration.transformOrder);
matlabInverseOutput = x;

for iRound = 1:(nWarmups+nSamples)
    firstOperation = mod(iRound-1,nOperations)+1;
    operationOrder = [firstOperation:nOperations 1:firstOperation-1];
    for iOperation = operationOrder
        if ~operationAvailable(iOperation)
            continue
        end
        try
            switch operationIds(iOperation)
                case "matlab-forward"
                    timer = tic;
                    matlabForwardOutput = fftAlong(x,configuration.transformOrder);
                    elapsed = toc(timer);
                    kernel = NaN;
                    pipeline = NaN;
                    inputAlignment = NaN;
                    outputAlignment = NaN;
                case "matlab-inverse"
                    timer = tic;
                    matlabInverseOutput = ifftAlong(matlabForwardOutput,configuration.transformOrder);
                    elapsed = toc(timer);
                    kernel = NaN;
                    pipeline = NaN;
                    inputAlignment = NaN;
                    outputAlignment = NaN;
                case "engine-forward-allocating"
                    timer = tic;
                    [allocatedOutput,kernel,pipeline,inputAlignment,outputAlignment] = feval(configuration.module,'r2c',plan,x);
                    elapsed = toc(timer);
                case "engine-forward-preallocated"
                    timer = tic;
                    [preallocatedSpectrum,kernel,pipeline,inputAlignment,outputAlignment] = feval(configuration.module,'r2c',plan,x,preallocatedSpectrum);
                    elapsed = toc(timer);
                case "engine-inverse-destructive"
                    destroyedSpectrum = complex(referenceSpectrum);
                    timer = tic;
                    [destroyedSpectrum,realOutput,kernel,pipeline,inputAlignment,outputAlignment] = feval(configuration.module,'c2r',plan,destroyedSpectrum,realOutput);
                    elapsed = toc(timer);
                    inverseOutput = realOutput;
            end
            if iRound > nWarmups
                sampleIndex = iRound-nWarmups;
                totalSamples(iOperation,sampleIndex) = elapsed;
                kernelSamples(iOperation,sampleIndex) = kernel;
                pipelineSamples(iOperation,sampleIndex) = pipeline;
                inputAlignmentSamples(iOperation,sampleIndex) = inputAlignment;
                outputAlignmentSamples(iOperation,sampleIndex) = outputAlignment;
            end
        catch exception
            operationAvailable(iOperation) = false;
            operationFailures(iOperation) = string(exception.identifier) + ": " + string(exception.message);
        end
    end
end

offset = 0;
measurement.matlabForwardSamplesSeconds = [];
measurement.matlabInverseSamplesSeconds = [];
measurement.matlabForwardMedianSeconds = NaN;
measurement.matlabInverseMedianSeconds = NaN;
measurement.matlabInverseRelativeError = NaN;
if includeMatlab
    measurement.matlabForwardSamplesSeconds = totalSamples(1,:);
    measurement.matlabInverseSamplesSeconds = totalSamples(2,:);
    measurement.matlabForwardMedianSeconds = median(totalSamples(1,:));
    measurement.matlabInverseMedianSeconds = median(totalSamples(2,:));
    [~,measurement.matlabInverseRelativeError] = numericalError(matlabInverseOutput,x);
    offset = 2;
else
    measurement.matlabInverseRelativeError = NaN;
end

measurement.forward.allocating = timingRecord(totalSamples(offset+1,:),kernelSamples(offset+1,:),pipelineSamples(offset+1,:),inputAlignmentSamples(offset+1,:),outputAlignmentSamples(offset+1,:),operationAvailable(offset+1),operationFailures(offset+1));
measurement.forward.preallocated = timingRecord(totalSamples(offset+2,:),kernelSamples(offset+2,:),pipelineSamples(offset+2,:),inputAlignmentSamples(offset+2,:),outputAlignmentSamples(offset+2,:),operationAvailable(offset+2),operationFailures(offset+2));
if measurement.forward.allocating.totalMexMedianSeconds <= measurement.forward.preallocated.totalMexMedianSeconds
    measurement.forward.winningAllocationMode = "allocating";
    winningForward = measurement.forward.allocating;
else
    measurement.forward.winningAllocationMode = "preallocated";
    winningForward = measurement.forward.preallocated;
end
measurement.forward.totalMexMedianSeconds = winningForward.totalMexMedianSeconds;
measurement.forward.kernelMedianSeconds = winningForward.kernelMedianSeconds;
measurement.forward.rawPipelineMedianSeconds = winningForward.rawPipelineMedianSeconds;
measurement.inverse = timingRecord(totalSamples(offset+3,:),kernelSamples(offset+3,:),pipelineSamples(offset+3,:),inputAlignmentSamples(offset+3,:),outputAlignmentSamples(offset+3,:),operationAvailable(offset+3),operationFailures(offset+3));

if operationAvailable(offset+1)
    [measurement.correctness.allocatingForwardAbsoluteError,measurement.correctness.allocatingForwardRelativeError] = numericalError(allocatedOutput,referenceSpectrum);
else
    measurement.correctness.allocatingForwardAbsoluteError = NaN;
    measurement.correctness.allocatingForwardRelativeError = NaN;
end
if operationAvailable(offset+2)
    [measurement.correctness.preallocatedForwardAbsoluteError,measurement.correctness.preallocatedForwardRelativeError] = numericalError(preallocatedSpectrum,referenceSpectrum);
else
    measurement.correctness.preallocatedForwardAbsoluteError = NaN;
    measurement.correctness.preallocatedForwardRelativeError = NaN;
end
scaleFactor = 1/prod(configuration.size(configuration.transformOrder));
if operationAvailable(offset+3)
    [measurement.correctness.inverseAbsoluteError,measurement.correctness.inverseRelativeError] = numericalError(scaleFactor*inverseOutput,x);
else
    measurement.correctness.inverseAbsoluteError = NaN;
    measurement.correctness.inverseRelativeError = NaN;
end
measurement.maximumRelativeError = max([measurement.correctness.allocatingForwardRelativeError measurement.correctness.preallocatedForwardRelativeError measurement.correctness.inverseRelativeError],[],'omitnan');
measurement.destructiveInputChanged = operationAvailable(offset+3) && ~isequaln(destroyedSpectrum,referenceSpectrum);
measurement.nWarmups = nWarmups;
measurement.nSamples = nSamples;
end

function record = timingRecord(totalSamples,kernelSamples,pipelineSamples,inputAlignmentSamples,outputAlignmentSamples,isAvailable,failure)
record.isAvailable = isAvailable;
record.failure = failure;
record.totalMexSamplesSeconds = totalSamples;
record.kernelSamplesSeconds = kernelSamples;
record.rawPipelineSamplesSeconds = pipelineSamples;
record.inputAlignmentClasses = inputAlignmentSamples;
record.outputAlignmentClasses = outputAlignmentSamples;
record.totalMexMedianSeconds = finiteMedian(totalSamples);
record.kernelMedianSeconds = finiteMedian(kernelSamples);
record.rawPipelineMedianSeconds = finiteMedian(pipelineSamples);
end

function value = finiteMedian(samples)
samples = samples(isfinite(samples));
if isempty(samples)
    value = Inf;
else
    value = median(samples);
end
end

function finalists = measureEngineFinalists(candidates,engine,options)
finalists = repmat(emptyWorkload(),0,1);
for iSize = 1:size(options.sizes,1)
    winner = selectWinner(candidates,[candidates.engine] == engine & [candidates.workloadIndex] == iSize);
    if isempty(winner)
        continue
    end
    finalists(end+1,1) = measureWinner(winner,options,iSize); %#ok<AGROW>
end
end

function workloads = measureGlobalFinalists(candidates,options)
workloads = repmat(emptyWorkload(),0,1);
for iSize = 1:size(options.sizes,1)
    winner = selectWinner(candidates,[candidates.workloadIndex] == iSize);
    if isempty(winner)
        error('FFTWEngineBenchmark:NoValidCandidate','No valid engine configuration remains for workload [%s].',num2str(options.sizes(iSize,:)));
    end
    workloads(end+1,1) = measureWinner(winner,options,iSize); %#ok<AGROW>
end
end

function winner = selectWinner(candidates,mask)
indices = find(mask & [candidates.status] == "passed");
if isempty(indices)
    winner = [];
    return
end
totalTimes = arrayfun(@(index) candidates(index).screening.forward.totalMexMedianSeconds,indices);
rawTimes = arrayfun(@(index) candidates(index).screening.forward.rawPipelineMedianSeconds,indices);
[~,order] = sortrows([totalTimes(:) rawTimes(:)],[1 2]);
winner = candidates(indices(order(1)));
end

function workload = measureWinner(winner,options,iSize)
rng(options.seed+iSize-1,'twister');
x = randn(winner.size);
[referenceSpectrum,complexSize] = spectrumReference(x,winner.transformOrder);
configuration = candidateConfiguration(winner.engine,winner.module,winner.size,winner.workloadIndex,winner.transformOrder,winner.strategy,winner.planner,winner.threads,winner.alignmentMode);
[plan,planRecord] = createPlan(configuration,x,complexSize,options.plannerTimeLimitSeconds);
planCleanup = onCleanup(@() feval(configuration.module,'free',plan));
nSamples = options.nSamples;
if iSize == size(options.sizes,1)
    nSamples = options.nSamplesLargest;
end
measurement = measureConfiguration(configuration,plan,x,referenceSpectrum,complexSize,options.nWarmups,nSamples,true);
alignedCeilingSamples = [];
if contains(winner.engine,"fftw")
    alignedCeilingSamples = feval(configuration.module,'alignedCeiling',plan,x,nSamples);
end
clear planCleanup

workload = emptyWorkload();
workload.size = winner.size;
workload.seed = options.seed+iSize-1;
workload.isGateWorkload = ismember(winner.size,options.gateSizes,'rows');
workload.winner = winner;
workload.winner.plan = planRecord;
workload.finalMeasurement = measurement;
workload.storage = storageRecord(winner.size,complexSize);
workload.rawForwardSpeedRatio = measurement.matlabForwardMedianSeconds/measurement.forward.rawPipelineMedianSeconds;
workload.totalForwardSpeedRatio = measurement.matlabForwardMedianSeconds/measurement.forward.totalMexMedianSeconds;
workload.inverseSpeedRatio = measurement.matlabInverseMedianSeconds/measurement.inverse.totalMexMedianSeconds;
workload.correctnessPassed = measurement.maximumRelativeError <= options.errorTolerance && measurement.destructiveInputChanged;
workload.rawForwardThresholdPassed = workload.rawForwardSpeedRatio >= options.rawForwardSpeedThreshold;
workload.totalForwardThresholdPassed = workload.totalForwardSpeedRatio >= options.totalForwardSpeedThreshold;
workload.inverseThresholdPassed = workload.inverseSpeedRatio >= options.inverseSpeedThreshold;
workload.fftwOwnedAlignedCeilingSamplesSeconds = alignedCeilingSamples;
if isempty(alignedCeilingSamples)
    workload.fftwOwnedAlignedCeilingMedianSeconds = NaN;
else
    workload.fftwOwnedAlignedCeilingMedianSeconds = median(alignedCeilingSamples);
end
end

function assessment = assessFinalists(finalists,options)
assessment = assessWorkloads(finalists,options);
end

function assessment = assessWorkloads(workloads,options)
assessment = emptyAssessment();
if isempty(workloads)
    return
end
gateMask = arrayfun(@(workload) ismember(workload.size,options.gateSizes,'rows'),workloads);
assessment.gateWorkloadCount = sum(gateMask);
assessment.requiredGateWorkloadCount = size(options.gateSizes,1);
assessment.wasEvaluated = assessment.gateWorkloadCount == assessment.requiredGateWorkloadCount;
if assessment.wasEvaluated
    assessment.rawForwardPassed = all([workloads(gateMask).rawForwardThresholdPassed]);
    assessment.totalForwardPassed = all([workloads(gateMask).totalForwardThresholdPassed]);
    assessment.forwardGatePassed = assessment.rawForwardPassed && assessment.totalForwardPassed;
    assessment.inversePassed = all([workloads(gateMask).inverseThresholdPassed]);
    assessment.correctnessPassed = all([workloads(gateMask).correctnessPassed]);
end
end

function [reference,complexSize] = spectrumReference(x,transformOrder)
fullSpectrum = fftAlong(x,transformOrder);
complexSize = size(x);
compressedDimension = transformOrder(2);
complexSize(compressedDimension) = floor(complexSize(compressedDimension)/2)+1;
indices = repmat({':'},1,ndims(x));
indices{compressedDimension} = 1:complexSize(compressedDimension);
reference = fullSpectrum(indices{:});
end

function output = fftAlong(input,transformOrder)
output = fft(input,[],transformOrder(1));
output = fft(output,[],transformOrder(2));
end

function output = ifftAlong(input,transformOrder)
output = ifft(input,[],transformOrder(2));
output = ifft(output,[],transformOrder(1),'symmetric');
end

function [maximumAbsoluteError,maximumRelativeError] = numericalError(actual,reference)
maximumAbsoluteError = max(abs(actual-reference),[],'all');
referenceMagnitude = max(abs(reference),[],'all');
if referenceMagnitude == 0, referenceMagnitude = 1; end
maximumRelativeError = maximumAbsoluteError/referenceMagnitude;
end

function storage = storageRecord(realSize,complexSize)
storage.realInputBytes = prod(realSize)*8;
storage.matlabFullSpectrumBytes = prod(realSize)*16;
storage.halfSpectrumBytes = prod(complexSize)*16;
storage.realOutputBytes = prod(realSize)*8;
storage.halfToFullSpectrumStorageRatio = storage.halfSpectrumBytes/storage.matlabFullSpectrumBytes;
end

function value = pointerAlignment(~)
value = "observed-by-native-gateway";
end

function value = layoutName(compressedDimension)
names = ["half-x","half-y","half-z"];
value = names(compressedDimension);
end

function environment = collectEnvironment(observedThreads)
environment.matlabVersion = string(version);
environment.matlabRelease = string(version('-release'));
environment.architecture = string(computer('arch'));
environment.mexExtension = string(mexext);
environment.operatingSystem = string(system_dependent('getos'));
environment.processor = "unknown";
environment.machineModel = "unknown";
environment.machineName = "unknown";
environment.physicalMemory = "unknown";
[status,text] = system('system_profiler SPHardwareDataType -json');
if status == 0
    data = jsondecode(text);
    hardware = data.SPHardwareDataType(1);
    environment.processor = optionalStructString(hardware,'chip_type',environment.processor);
    environment.machineModel = optionalStructString(hardware,'machine_model',environment.machineModel);
    environment.machineName = optionalStructString(hardware,'machine_name',environment.machineName);
    environment.physicalMemory = optionalStructString(hardware,'physical_memory',environment.physicalMemory);
end
environment.observedComputationalThreads = observedThreads;
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
if isfield(inputStruct,fieldName), value = string(inputStruct.(fieldName)); else, value = defaultValue; end
end

function sources = collectSourceProvenance(sourceDirectory)
names = ["runFFTWEngineLayoutBenchmark.m","buildFFTWEngineBenchmarkMex.m","fftw_engine_benchmark.cpp","vdsp_engine_benchmark.cpp"];
files = repmat(struct('path',"",'sha256',""),numel(names),1);
for iFile = 1:numel(names)
    files(iFile).path = names(iFile);
    files(iFile).sha256 = textFileSHA256(fullfile(sourceDirectory,names(iFile)));
end
sources.files = files;
repositoryRoot = fileparts(fileparts(fileparts(sourceDirectory)));
[status,commit] = system(sprintf('git -C "%s" rev-parse HEAD',repositoryRoot));
if status == 0, sources.repositoryCommit = string(strtrim(commit)); else, sources.repositoryCommit = "unknown"; end
[status,dirty] = system(sprintf('git -C "%s" status --porcelain --untracked-files=all',repositoryRoot));
sources.repositoryDirty = status ~= 0 || strlength(strtrim(string(dirty))) > 0;
end

function hash = textFileSHA256(path)
digest = java.security.MessageDigest.getInstance('SHA-256');
digest.update(unicode2native(fileread(path),'UTF-8'));
hashBytes = typecast(digest.digest(),'uint8');
hash = string(lower(reshape(dec2hex(hashBytes,2).',1,[])));
end

function writeResultArtifacts(result,runDirectory)
serializable = result;
serializable.engines = num2cell(result.engines);
serializable.alignmentValidation = num2cell(result.alignmentValidation);
serializable.candidates = num2cell(result.candidates);
serializable.workloads = num2cell(result.workloads);
jsonText = jsonencode(serializable);
writeTextFile(fullfile(runDirectory,"engine-layout-benchmark.json"),string(jsonText));
writeTextFile(fullfile(runDirectory,"summary.md"),createMarkdownSummary(result));
end

function summary = createMarkdownSummary(result)
lines = strings(0,1);
lines(end+1) = "# FFT engine and half-spectrum layout benchmark";
lines(end+1) = "";
lines(end+1) = "Status: **" + upper(result.status) + "**";
lines(end+1) = "";
lines(end+1) = "This issue #38 report shows criterion results but does not make the milestone's GO or NO-GO decision.";
lines(end+1) = "";
lines(end+1) = "## Environment";
lines(end+1) = "";
lines(end+1) = "| Field | Value |";
lines(end+1) = "|---|---|";
lines(end+1) = markdownRow("Run",result.runId);
lines(end+1) = markdownRow("MATLAB",result.environment.matlabVersion);
lines(end+1) = markdownRow("Architecture",result.environment.architecture);
lines(end+1) = markdownRow("Operating system",result.environment.operatingSystem);
lines(end+1) = markdownRow("Processor",result.environment.processor);
lines(end+1) = markdownRow("Machine",result.environment.machineName + " " + result.environment.machineModel);
lines(end+1) = markdownRow("Physical memory",result.environment.physicalMemory);
lines(end+1) = markdownRow("Threads",string(result.environment.observedComputationalThreads));
lines(end+1) = markdownRow("Planner limit",string(result.configuration.plannerTimeLimitSeconds) + " seconds per plan");

lines(end+1) = "";
lines(end+1) = "## Engines";
lines(end+1) = "";
lines(end+1) = "| Engine | Status | Version | Library | Reason |";
lines(end+1) = "|---|---|---|---|---|";
for engine = result.engines'
    lines(end+1) = "| " + escapeMarkdown(engine.id) + " | " + engine.status + " | " + escapeMarkdown(engine.version) + " | " + escapeMarkdown(engine.library) + " | " + escapeMarkdown(engine.reason) + " |"; %#ok<AGROW>
end
if ~isempty(result.alignmentValidation)
    lines(end+1) = "";
    lines(end+1) = "## FFTW new-array alignment validation";
    lines(end+1) = "";
    lines(end+1) = "| Engine | Matched accepted | Mismatch rejected | Distinct offset classes observed |";
    lines(end+1) = "|---|:---:|:---:|:---:|";
    for validation = result.alignmentValidation'
        lines(end+1) = sprintf('| %s | %s | %s | %s |',validation.engine,yesNo(validation.matchedAccepted),yesNo(validation.mismatchRejected),yesNo(validation.distinctPointerClassesObserved)); %#ok<AGROW>
    end
end

if ~isempty(result.workloads)
    lines(end+1) = "";
    lines(end+1) = "## Fastest valid configuration by workload";
    lines(end+1) = "";
    lines(end+1) = "Speed ratios above 1 are faster than the corresponding complete MATLAB transform.";
    lines(end+1) = "";
    lines(end+1) = "| Size | Engine | Layout | Strategy | Planner | Threads | Alignment | Output | Raw r2c | Total MEX r2c | Destructive c2r | Relative error | Half/full |";
    lines(end+1) = "|---|---|---|---|---|---:|---|---|---:|---:|---:|---:|---:|";
    for workload = result.workloads'
        winner = workload.winner;
        measurement = workload.finalMeasurement;
        lines(end+1) = sprintf('| %s | %s | %s | %s | %s | %d | %s | %s | %.3fx | %.3fx | %.3fx | %.3g | %.6f |',formatSize(workload.size),winner.engine,winner.layout,winner.strategy,winner.planner,winner.threads,winner.alignmentMode,measurement.forward.winningAllocationMode,workload.rawForwardSpeedRatio,workload.totalForwardSpeedRatio,workload.inverseSpeedRatio,measurement.maximumRelativeError,workload.storage.halfToFullSpectrumStorageRatio); %#ok<AGROW>
    end

    lines(end+1) = "";
    lines(end+1) = "## Criterion status";
    lines(end+1) = "";
    lines(end+1) = "The performance gate applies only to `256 x 256 x 64` and `1024 x 1024 x 30`.";
    lines(end+1) = "";
    lines(end+1) = "| Size | Gate workload | Raw r2c >= 1.25x | Total MEX r2c >= 1.10x | Destructive c2r >= 0.95x | Error <= 1e-12 |";
    lines(end+1) = "|---|:---:|:---:|:---:|:---:|:---:|";
    for workload = result.workloads'
        lines(end+1) = sprintf('| %s | %s | %s | %s | %s | %s |',formatSize(workload.size),yesNo(workload.isGateWorkload),yesNo(workload.rawForwardThresholdPassed),yesNo(workload.totalForwardThresholdPassed),yesNo(workload.inverseThresholdPassed),yesNo(workload.correctnessPassed)); %#ok<AGROW>
    end
end
lines(end+1) = "";
lines(end+1) = "## Timing boundaries";
lines(end+1) = "";
lines(end+1) = "- Kernel time contains only FFTW execute calls or vDSP transform calls.";
lines(end+1) = "- Raw pipeline time also contains required vDSP packing, unpacking, scaling, and canonical layout conversion.";
lines(end+1) = "- Total MEX time contains the complete MATLAB-to-MEX call using the existing allocation models. Detailed ownership attribution remains issue #39.";
if result.status == "failed" && ~isempty(result.failure)
    lines(end+1) = "";
    lines(end+1) = "## Failure";
    lines(end+1) = "";
    lines(end+1) = "- Stage: `" + result.failure.stage + "`";
    lines(end+1) = "- Identifier: `" + result.failure.identifier + "`";
    lines(end+1) = "- Message: " + escapeMarkdown(result.failure.message);
end
summary = strjoin(lines,newline);
end

function row = markdownRow(field,value)
row = "| " + escapeMarkdown(string(field)) + " | " + escapeMarkdown(string(value)) + " |";
end

function value = escapeMarkdown(value)
value = replace(string(value),"|","\|");
value = replace(value,newline,"<br>");
end

function value = formatSize(sz)
value = strtrim(sprintf('%d x %d x %d',sz));
end

function value = yesNo(tf)
if tf, value = 'yes'; else, value = 'no'; end
end

function writeTextFile(path,text)
[fileId,message] = fopen(path,'w');
if fileId < 0, error('FFTWEngineBenchmark:ArtifactOpenFailed','Unable to open %s: %s',path,message); end
cleanup = onCleanup(@() fclose(fileId));
fprintf(fileId,'%s\n',text);
clear cleanup
end

function restoreMatlabState(previousPlanner,previousWisdom,previousThreads,previousRandomState)
try
    fftw('dwisdom',[]); fftw('dwisdom',previousWisdom);
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

function record = exceptionRecord(exception)
record.identifier = string(exception.identifier);
record.message = string(exception.message);
end

function failure = failureRecord(exception,stage)
failure = exceptionRecord(exception);
failure.stage = stage;
failure.stack = repmat(struct('name',"",'line',0),numel(exception.stack),1);
for iFrame = 1:numel(exception.stack)
    failure.stack(iFrame).name = string(exception.stack(iFrame).name);
    failure.stack(iFrame).line = exception.stack(iFrame).line;
end
end

function timestamp = utcTimestamp()
timestamp = string(datetime('now',TimeZone='UTC',Format="yyyy-MM-dd'T'HH:mm:ss.SSS'Z'"));
end

function engine = availableEngine(id,module,version,library)
engine = emptyEngine();
engine.id = id;
engine.module = module;
engine.status = "available";
engine.version = version;
engine.library = library;
end

function engine = unavailableEngine(id,exception)
engine = unavailableEngineMessage(id,string(exception.identifier) + ": " + string(exception.message));
end

function engine = unavailableEngineMessage(id,message)
engine = emptyEngine();
engine.id = id;
engine.status = "unavailable";
engine.reason = message;
end

function engine = emptyEngine()
engine.id = "";
engine.module = "";
engine.status = "";
engine.version = "";
engine.library = "";
engine.reason = "";
end

function validation = alignmentValidationRecord(engine,module)
[matchedAccepted,mismatchRejected,distinctPointerClassesObserved] = feval(module,'alignmentSelfTest');
validation = emptyAlignmentValidation();
validation.engine = engine;
validation.matchedAccepted = logical(matchedAccepted);
validation.mismatchRejected = logical(mismatchRejected);
validation.distinctPointerClassesObserved = logical(distinctPointerClassesObserved);
end

function validation = emptyAlignmentValidation()
validation.engine = "";
validation.matchedAccepted = false;
validation.mismatchRejected = false;
validation.distinctPointerClassesObserved = false;
end

function candidate = emptyCandidate()
candidate.engine = "";
candidate.module = "";
candidate.size = [];
candidate.workloadIndex = 0;
candidate.transformOrder = [];
candidate.layout = "";
candidate.strategy = "";
candidate.planner = "";
candidate.threads = 0;
candidate.alignmentMode = "";
candidate.status = "";
candidate.plan = struct('planningSeconds',NaN,'planningLimitReached',false,'inputAlignmentClass',NaN,'outputAlignmentClass',NaN,'observedThreads',"");
candidate.screening = emptyMeasurement();
candidate.storage = struct('realInputBytes',0,'matlabFullSpectrumBytes',0,'halfSpectrumBytes',0,'realOutputBytes',0,'halfToFullSpectrumStorageRatio',NaN);
candidate.maximumRelativeError = NaN;
candidate.destructiveInputChanged = false;
candidate.failure = struct('identifier',"",'message',"");
end

function measurement = emptyMeasurement()
emptyTiming = timingRecord([],[],[],[],[],false,"");
measurement.matlabForwardSamplesSeconds = [];
measurement.matlabInverseSamplesSeconds = [];
measurement.matlabForwardMedianSeconds = NaN;
measurement.matlabInverseMedianSeconds = NaN;
measurement.forward.allocating = emptyTiming;
measurement.forward.preallocated = emptyTiming;
measurement.forward.winningAllocationMode = "";
measurement.forward.totalMexMedianSeconds = Inf;
measurement.forward.kernelMedianSeconds = Inf;
measurement.forward.rawPipelineMedianSeconds = Inf;
measurement.inverse = emptyTiming;
measurement.correctness = struct('allocatingForwardAbsoluteError',NaN,'allocatingForwardRelativeError',NaN,'preallocatedForwardAbsoluteError',NaN,'preallocatedForwardRelativeError',NaN,'inverseAbsoluteError',NaN,'inverseRelativeError',NaN);
measurement.maximumRelativeError = NaN;
measurement.destructiveInputChanged = false;
measurement.nWarmups = 0;
measurement.nSamples = 0;
end

function workload = emptyWorkload()
workload.size = [];
workload.seed = 0;
workload.isGateWorkload = false;
workload.winner = emptyCandidate();
workload.finalMeasurement = emptyMeasurement();
workload.storage = struct('realInputBytes',0,'matlabFullSpectrumBytes',0,'halfSpectrumBytes',0,'realOutputBytes',0,'halfToFullSpectrumStorageRatio',NaN);
workload.rawForwardSpeedRatio = NaN;
workload.totalForwardSpeedRatio = NaN;
workload.inverseSpeedRatio = NaN;
workload.correctnessPassed = false;
workload.rawForwardThresholdPassed = false;
workload.totalForwardThresholdPassed = false;
workload.inverseThresholdPassed = false;
workload.fftwOwnedAlignedCeilingSamplesSeconds = [];
workload.fftwOwnedAlignedCeilingMedianSeconds = NaN;
end

function assessment = emptyAssessment()
assessment.wasEvaluated = false;
assessment.gateWorkloadCount = 0;
assessment.requiredGateWorkloadCount = 0;
assessment.rawForwardPassed = false;
assessment.totalForwardPassed = false;
assessment.forwardGatePassed = false;
assessment.inversePassed = false;
assessment.correctnessPassed = false;
end
