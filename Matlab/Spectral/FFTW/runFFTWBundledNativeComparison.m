function result = runFFTWBundledNativeComparison(options)
% Compare MATLAB-bundled and native FFTW with identical configurations.
%
% This issue #41 benchmark determines horizontal bundled-FFTW readiness. It
% uses native FFTW only as a controlled reference and does not decide the
% readiness of real-to-real transforms or the WaveVortex adapter.
arguments
    options.outputDirectory (1,1) string = fullfile(fileparts(mfilename('fullpath')),"results","issue41")
    options.issue38ArtifactPath (1,1) string = fullfile(fileparts(mfilename('fullpath')),"results","issue38","20260808T045405991Z-maca64-r2026a","engine-layout-benchmark.json")
    options.sizes (:,3) double = [256 256 64; 1024 1024 30]
    options.gateSizes (:,3) double = [256 256 64; 1024 1024 30]
    options.transformOrder (1,2) double = [2 1]
    options.planners (1,:) string = ["estimate","measure","patient","exhaustive"]
    options.alignmentModes (1,:) string = ["matched","unaligned"]
    options.threadCount (1,1) double = maxNumCompThreads
    options.plannerTimeLimitSeconds (1,1) double = 10
    options.nWarmups (1,1) double = 2
    options.nSamples (1,1) double = 7
    options.nSamplesLargest (1,1) double = 3
    options.seed (1,1) double = 41
    options.errorTolerance (1,1) double = 1e-12
    options.rawForwardSpeedThreshold (1,1) double = 1.25
    options.totalForwardSpeedThreshold (1,1) double = 1.10
    options.inverseSpeedThreshold (1,1) double = 0.95
    options.shouldReplayIssue38 (1,1) logical = true
    options.shouldBuild (1,1) logical = true
    options.shouldWriteArtifacts (1,1) logical = true
    options.runId (1,1) string = ""
end

validateOptions(options);
sourceDirectory = fileparts(mfilename('fullpath'));
runId = options.runId;
if strlength(runId) == 0
    runId = string(datetime('now',TimeZone='UTC',Format="yyyyMMdd'T'HHmmssSSS'Z'")) + "-" + string(computer('arch')) + "-r" + lower(string(version('-release')));
end
runDirectory = prepareRunDirectory(options.outputDirectory,runId,options.shouldWriteArtifacts);

previousPlanner = string(fftw('planner'));
previousWisdom = fftw('dwisdom');
previousThreads = maxNumCompThreads;
previousRandomState = rng;
stateCleanup = onCleanup(@() restoreMatlabState(previousPlanner,previousWisdom,previousThreads,previousRandomState));
maxNumCompThreads('automatic');
fftw('planner','measure');
fftw('dwisdom',[]);

result = initializeResult(runId,runDirectory,sourceDirectory,options);
activeStage = "build";
try
    if options.shouldBuild
        result.build = buildFFTWBundledNativeComparisonMex;
    else
        requireModules(options.shouldReplayIssue38);
        result.build = existingBuildRecord;
    end
    result.engines = collectEngineRecords;

    activeStage = "matched bundled/native measurements";
    requestedThreads = min(options.threadCount,maxNumCompThreads);
    result.configuration.requestedThreads = options.threadCount;
    result.configuration.observedThreads = requestedThreads;
    for iSize = 1:size(options.sizes,1)
        nSamples = options.nSamples;
        if iSize == size(options.sizes,1), nSamples = options.nSamplesLargest; end
        rng(options.seed+iSize-1,'twister');
        x = randn(options.sizes(iSize,:));
        for planner = options.planners
            for alignmentMode = options.alignmentModes
                configuration = configurationRecord(options.sizes(iSize,:),options.transformOrder,planner,requestedThreads,alignmentMode);
                try
                    comparison = benchmarkMatchedConfiguration(configuration,x,options.nWarmups,nSamples,options.plannerTimeLimitSeconds,options.errorTolerance);
                catch exception
                    comparison = failedComparison(configuration,options.nWarmups,nSamples,exception);
                end
                result.comparisons(end+1,1) = comparison; %#ok<AGROW>
            end
        end
    end

    activeStage = "workload readiness selection";
    result.workloads = selectWorkloadResults(result.comparisons,options);
    result.readiness = readinessRecord(result.workloads,options);

    if options.shouldReplayIssue38
        activeStage = "issue #38 historical replay";
        result.historicalReplay = replayIssue38(options,result.workloads);
    end

    result.status = "passed";
    result.completedAtUTC = utcTimestamp;
    if options.shouldWriteArtifacts, writeResultArtifacts(result,runDirectory); end
catch exception
    result.status = "failed";
    result.completedAtUTC = utcTimestamp;
    result.failure = failureRecord(exception,activeStage);
    if options.shouldWriteArtifacts
        try
            writeResultArtifacts(result,runDirectory);
        catch artifactException
            warning('FFTWBundledNative:FailureArtifactWriteFailed','Unable to write issue #41 failure artifacts: %s',artifactException.message);
        end
    end
    rethrow(exception);
end
clear stateCleanup
end

function validateOptions(options)
validateattributes(options.sizes,{'double'},{'2d','ncols',3,'integer','positive'},mfilename,'sizes');
validateattributes(options.gateSizes,{'double'},{'2d','ncols',3,'integer','positive'},mfilename,'gateSizes');
validateattributes(options.transformOrder,{'double'},{'numel',2,'integer','positive'},mfilename,'transformOrder');
if options.transformOrder(1) == options.transformOrder(2) || any(options.transformOrder > 3)
    error('FFTWBundledNative:InvalidTransformOrder','The transform order must contain two distinct dimensions between 1 and 3.');
end
if any(~ismember(options.planners,["estimate","measure","patient","exhaustive"]))
    error('FFTWBundledNative:UnknownPlanner','Planners must be estimate, measure, patient, or exhaustive.');
end
if any(~ismember(options.alignmentModes,["matched","unaligned"]))
    error('FFTWBundledNative:UnknownAlignmentMode','Alignment modes must be matched or unaligned.');
end
validateattributes(options.threadCount,{'double'},{'scalar','integer','positive'},mfilename,'threadCount');
validateattributes(options.plannerTimeLimitSeconds,{'double'},{'scalar','positive','finite'},mfilename,'plannerTimeLimitSeconds');
validateattributes(options.nWarmups,{'double'},{'scalar','integer','nonnegative'},mfilename,'nWarmups');
validateattributes(options.nSamples,{'double'},{'scalar','integer','positive'},mfilename,'nSamples');
validateattributes(options.nSamplesLargest,{'double'},{'scalar','integer','positive'},mfilename,'nSamplesLargest');
validateattributes(options.seed,{'double'},{'scalar','integer','nonnegative'},mfilename,'seed');
validateattributes(options.errorTolerance,{'double'},{'scalar','real','finite','nonnegative'},mfilename,'errorTolerance');
if ~startsWith(string(version('-release')),"2026a",IgnoreCase=true) || string(computer('arch')) ~= "maca64"
    error('FFTWBundledNative:UnsupportedPlatform','Issue #41 targets MATLAB R2026a on macOS maca64.');
end
end

function runDirectory = prepareRunDirectory(outputDirectory,runId,shouldWrite)
runDirectory = "";
if ~shouldWrite, return; end
runDirectory = fullfile(outputDirectory,runId);
if isfolder(runDirectory) || isfile(runDirectory)
    error('FFTWBundledNative:OutputExists','Benchmark output already exists at %s.',runDirectory);
end
[didCreate,message] = mkdir(runDirectory);
if ~didCreate, error('FFTWBundledNative:OutputCreationFailed','Unable to create %s: %s',runDirectory,message); end
end

function result = initializeResult(runId,runDirectory,sourceDirectory,options)
result.schemaVersion = "1.0.0";
result.status = "running";
result.runId = runId;
result.generatedAtUTC = utcTimestamp;
result.completedAtUTC = "";
result.environment = collectEnvironment;
result.configuration.sizes = options.sizes;
result.configuration.gateSizes = options.gateSizes;
result.configuration.transformOrder = options.transformOrder;
result.configuration.layout = layoutName(options.transformOrder(2));
result.configuration.strategy = "guru-rank2";
result.configuration.planners = options.planners;
result.configuration.alignmentModes = options.alignmentModes;
result.configuration.requestedThreads = options.threadCount;
result.configuration.observedThreads = NaN;
result.configuration.plannerTimeLimitSeconds = options.plannerTimeLimitSeconds;
result.configuration.nWarmups = options.nWarmups;
result.configuration.nSamples = options.nSamples;
result.configuration.nSamplesLargest = options.nSamplesLargest;
result.configuration.seeds = options.seed + (0:size(options.sizes,1)-1);
result.configuration.errorTolerance = options.errorTolerance;
result.configuration.thresholds.rawForwardSpeedRatio = options.rawForwardSpeedThreshold;
result.configuration.thresholds.totalForwardSpeedRatio = options.totalForwardSpeedThreshold;
result.configuration.thresholds.inverseSpeedRatio = options.inverseSpeedThreshold;
result.configuration.forwardOwnership = "MATLAB createBuffer/createArrayFromBuffer";
result.configuration.inverseOwnership = "uniquely owned destructive input and caller-preallocated real output";
result.configuration.selection = "all matched pairs fully sampled; lowest bundled total forward MEX median, raw kernel tiebreaker";
result.configuration.operationOrder = "six-operation round-robin with rotating first operation";
result.configuration.issue38ArtifactPath = options.issue38ArtifactPath;
result.sources = collectSourceProvenance(sourceDirectory);
result.build = struct;
result.engines = repmat(emptyEngine,0,1);
result.comparisons = repmat(emptyComparison,0,1);
result.workloads = repmat(emptyWorkload,0,1);
result.historicalReplay = repmat(emptyHistoricalReplay,0,1);
result.readiness = emptyReadiness;
result.failure = [];
result.artifacts.directory = runDirectory;
result.artifacts.json = "bundled-native-comparison.json";
result.artifacts.markdown = "summary.md";
end

function requireModules(includeReplay)
modules = ["fftw_ownership_benchmark_bundled","fftw_ownership_benchmark_native"];
if includeReplay, modules = [modules,"fftw_engine_benchmark_bundled","fftw_engine_benchmark_native"]; end
for module = modules
    if exist(module,'file') ~= 3
        error('FFTWBundledNative:MexMissing','%s is unavailable and shouldBuild is false.',module);
    end
end
end

function build = existingBuildRecord()
build.status = "existing";
build.ownership.bundled.module = "fftw_ownership_benchmark_bundled";
build.ownership.native.module = "fftw_ownership_benchmark_native";
build.replay.bundled.module = "fftw_engine_benchmark_bundled";
build.replay.native.module = "fftw_engine_benchmark_native";
end

function engines = collectEngineRecords()
ids = ["bundled-fftw","native-fftw"];
modules = ["fftw_ownership_benchmark_bundled","fftw_ownership_benchmark_native"];
engines = repmat(emptyEngine,numel(ids),1);
for iEngine = 1:numel(ids)
    [engineVersion,library] = feval(modules(iEngine),'info');
    [matchedAccepted,mismatchRejected,unalignedAccepted] = feval(modules(iEngine),'alignmentSelfTest');
    engines(iEngine).id = ids(iEngine);
    engines(iEngine).module = modules(iEngine);
    engines(iEngine).version = string(engineVersion);
    engines(iEngine).library = string(library);
    engines(iEngine).alignmentMatchedAccepted = logical(matchedAccepted);
    engines(iEngine).alignmentMismatchRejected = logical(mismatchRejected);
    engines(iEngine).unalignedAccepted = logical(unalignedAccepted);
end
if startsWith(engines(2).library,string(matlabroot)) || engines(1).library == engines(2).library
    error('FFTWBundledNative:NativeResolutionFailure','The native comparison module resolved FFTW symbols to MATLAB''s bundled library.');
end
end

function configuration = configurationRecord(sz,transformOrder,planner,threads,alignmentMode)
configuration.size = sz;
configuration.transformOrder = transformOrder;
configuration.layout = layoutName(transformOrder(2));
configuration.strategy = "guru-rank2";
configuration.planner = planner;
configuration.threads = threads;
configuration.alignmentMode = alignmentMode;
configuration.forwardOwnership = "matlab-buffer";
configuration.inverseOwnership = "destructive-unique";
end

function comparison = benchmarkMatchedConfiguration(configuration,x,nWarmups,nSamples,timeLimit,errorTolerance)
modules = ["fftw_ownership_benchmark_bundled","fftw_ownership_benchmark_native"];
for module = modules, feval(module,'forgetWisdom'); end

[referenceSpectrum,complexSize] = spectrumReference(x,configuration.transformOrder);
realTemplate = zeros(configuration.size);
complexTemplate = complex(zeros(complexSize));
flags = plannerFlags(configuration.planner);
plans = zeros(1,2,'uint64');
planRecords = repmat(emptyPlan,2,1);
planCleanup = onCleanup(@() freePlans(modules,plans));
for iEngine = 1:2
    [plans(iEngine),reportedSize,scaleFactor,planningSeconds,inputAlignment,outputAlignment,limitReached] = feval(modules(iEngine),'create',configuration.size,configuration.transformOrder,configuration.threads,flags,char(configuration.alignmentMode),realTemplate,complexTemplate,timeLimit);
    if ~isequal(double(reportedSize),double(complexSize))
        error('FFTWBundledNative:ReportedSizeMismatch','%s reported an unexpected half-spectrum shape.',modules(iEngine));
    end
    planRecords(iEngine).planningSeconds = planningSeconds;
    planRecords(iEngine).planningLimitReached = logical(limitReached);
    planRecords(iEngine).inputAlignmentClass = inputAlignment;
    planRecords(iEngine).outputAlignmentClass = outputAlignment;
    planRecords(iEngine).wisdomClearedBeforePlanning = true;
end

operationIds = ["matlab-forward","matlab-inverse","bundled-forward","bundled-inverse","native-forward","native-inverse"];
nOperations = numel(operationIds);
totalSamples = nan(nOperations,nSamples);
metricSamples = nan(nOperations,17,nSamples);
pointerSamples = zeros(nOperations,5,nSamples,'uint64');
returnedTokens = zeros(nOperations,nSamples,'uint64');
matlabSpectrum = fftAlong(x,configuration.transformOrder);
matlabOutput = x;
bundledForward = complexTemplate;
nativeForward = complexTemplate;
bundledInverse = realTemplate;
nativeInverse = realTemplate;
destroyedBundled = referenceSpectrum;
destroyedNative = referenceSpectrum;

for iRound = 1:(nWarmups+nSamples)
    firstOperation = mod(iRound-1,nOperations)+1;
    operationOrder = [firstOperation:nOperations 1:firstOperation-1];
    for iOperation = operationOrder
        metrics = nan(1,17);
        pointers = zeros(1,5,'uint64');
        returned = uint64(0);
        switch operationIds(iOperation)
            case "matlab-forward"
                timer = tic;
                matlabSpectrum = fftAlong(x,configuration.transformOrder);
                elapsed = toc(timer);
            case "matlab-inverse"
                timer = tic;
                matlabOutput = ifftAlong(matlabSpectrum,configuration.transformOrder);
                elapsed = toc(timer);
            case "bundled-forward"
                timer = tic;
                bundledForward = feval(modules(1),'forward',plans(1),x,'matlab-buffer');
                elapsed = toc(timer);
                [metrics,pointers] = feval(modules(1),'metrics',plans(1));
                returned = feval(modules(1),'pointer',bundledForward);
            case "bundled-inverse"
                destroyedBundled = complex(real(referenceSpectrum),imag(referenceSpectrum));
                timer = tic;
                [destroyedBundled,bundledInverse] = feval(modules(1),'inverseDestructive',plans(1),destroyedBundled,bundledInverse);
                elapsed = toc(timer);
                [metrics,pointers] = feval(modules(1),'metrics',plans(1));
                returned = feval(modules(1),'pointer',destroyedBundled);
            case "native-forward"
                timer = tic;
                nativeForward = feval(modules(2),'forward',plans(2),x,'matlab-buffer');
                elapsed = toc(timer);
                [metrics,pointers] = feval(modules(2),'metrics',plans(2));
                returned = feval(modules(2),'pointer',nativeForward);
            case "native-inverse"
                destroyedNative = complex(real(referenceSpectrum),imag(referenceSpectrum));
                timer = tic;
                [destroyedNative,nativeInverse] = feval(modules(2),'inverseDestructive',plans(2),destroyedNative,nativeInverse);
                elapsed = toc(timer);
                [metrics,pointers] = feval(modules(2),'metrics',plans(2));
                returned = feval(modules(2),'pointer',destroyedNative);
        end
        if iRound > nWarmups
            iSample = iRound-nWarmups;
            totalSamples(iOperation,iSample) = elapsed;
            metricSamples(iOperation,:,iSample) = metrics;
            pointerSamples(iOperation,:,iSample) = pointers;
            returnedTokens(iOperation,iSample) = returned;
        end
    end
end

comparison = emptyComparison;
comparison.status = "passed";
comparison.configuration = configuration;
comparison.nWarmups = nWarmups;
comparison.nSamples = nSamples;
comparison.matlab.forward = matlabTiming(totalSamples(1,:));
comparison.matlab.inverse = matlabTiming(totalSamples(2,:));
comparison.matlab.maximumRelativeError = numericalError(matlabOutput,x);
for iEngine = 1:2
    forwardIndex = 2*iEngine+1;
    inverseIndex = forwardIndex+1;
    if iEngine == 1
        forwardOutput = bundledForward;
        inverseOutput = bundledInverse;
    else
        forwardOutput = nativeForward;
        inverseOutput = nativeInverse;
    end
    [forwardAbsolute,forwardRelative] = numericalErrors(forwardOutput,referenceSpectrum);
    [inverseAbsolute,inverseRelative] = numericalErrors(scaleFactor*inverseOutput,x);
    forward = engineTiming(totalSamples(forwardIndex,:),reshape(metricSamples(forwardIndex,:,:),17,nSamples).',reshape(pointerSamples(forwardIndex,:,:),5,nSamples).',returnedTokens(forwardIndex,:),forwardAbsolute,forwardRelative,errorTolerance,"forward");
    inverse = engineTiming(totalSamples(inverseIndex,:),reshape(metricSamples(inverseIndex,:,:),17,nSamples).',reshape(pointerSamples(inverseIndex,:,:),5,nSamples).',returnedTokens(inverseIndex,:),inverseAbsolute,inverseRelative,errorTolerance,"inverse");
    engine = emptyEngineMeasurement;
    engine.plan = planRecords(iEngine);
    engine.forward = forward;
    engine.inverse = inverse;
    engine.maximumRelativeError = max(forward.maximumRelativeError,inverse.maximumRelativeError);
    engine.correctnessPassed = forward.correctnessPassed && inverse.correctnessPassed;
    engine.zeroCopyPassed = forward.zeroCopyPassed && inverse.zeroCopyPassed;
    engine.destructiveInputChanged = all(inverse.destroyedInputSamples);
    engine.rawForwardSpeedRatio = comparison.matlab.forward.medianSeconds/forward.kernelMedianSeconds;
    engine.totalForwardSpeedRatio = comparison.matlab.forward.medianSeconds/forward.totalMexMedianSeconds;
    engine.inverseSpeedRatio = comparison.matlab.inverse.medianSeconds/inverse.totalMexMedianSeconds;
    if iEngine == 1, comparison.bundled = engine; else, comparison.native = engine; end
end
comparison.bundledRelativeToNative.forwardKernelPercent = percentDifference(comparison.bundled.forward.kernelMedianSeconds,comparison.native.forward.kernelMedianSeconds);
comparison.bundledRelativeToNative.forwardTotalPercent = percentDifference(comparison.bundled.forward.totalMexMedianSeconds,comparison.native.forward.totalMexMedianSeconds);
comparison.bundledRelativeToNative.inverseTotalPercent = percentDifference(comparison.bundled.inverse.totalMexMedianSeconds,comparison.native.inverse.totalMexMedianSeconds);
comparison.storage = storageRecord(configuration.size,complexSize);
clear planCleanup
end

function timing = matlabTiming(samples)
timing.samplesSeconds = samples;
timing.medianSeconds = median(samples);
end

function timing = engineTiming(totalSamples,metricSamples,pointerSamples,returnedTokens,absoluteError,relativeError,errorTolerance,direction)
timing = emptyEngineTiming;
timing.totalMexSamplesSeconds = totalSamples;
timing.totalMexMedianSeconds = median(totalSamples);
names = ["allocation","bufferWrap","explicitMemcpy","mutableDetach","kernel","internalPipeline"];
for iName = 1:numel(names)
    timing.(names(iName)+"SamplesSeconds") = metricSamples(:,iName).';
    timing.(names(iName)+"MedianSeconds") = median(metricSamples(:,iName));
end
timing.boundaryResidualSamplesSeconds = max(0,totalSamples-metricSamples(:,6).');
timing.boundaryResidualMedianSeconds = median(timing.boundaryResidualSamplesSeconds);
timing.explicitCopyCountSamples = metricSamples(:,9).';
timing.explicitCopiedBytesSamples = metricSamples(:,10).';
timing.detectedCopyCountSamples = metricSamples(:,11).';
timing.detectedCopiedBytesSamples = metricSamples(:,12).';
timing.inputAlignmentClasses = metricSamples(:,13).';
timing.outputAlignmentClasses = metricSamples(:,14).';
timing.destroyedInputSamples = logical(metricSamples(:,15).');
timing.inputBeforeTokens = pointerSamples(:,1).';
timing.inputMutableTokens = pointerSamples(:,2).';
timing.outputMutableTokens = pointerSamples(:,4).';
timing.wrappedTokens = pointerSamples(:,5).';
timing.returnedTokens = returnedTokens;
timing.maximumAbsoluteError = absoluteError;
timing.maximumRelativeError = relativeError;
timing.correctnessPassed = relativeError <= errorTolerance;
if direction == "forward"
    timing.pointerPreserved = all(timing.returnedTokens == timing.wrappedTokens);
    timing.zeroCopyPassed = timing.pointerPreserved && all(timing.explicitCopyCountSamples == 0) && all(timing.detectedCopyCountSamples == 0);
else
    timing.pointerPreserved = all(timing.returnedTokens == timing.inputMutableTokens);
    timing.zeroCopyPassed = timing.pointerPreserved && all(timing.inputBeforeTokens == timing.inputMutableTokens) && all(timing.explicitCopyCountSamples == 0) && all(timing.detectedCopyCountSamples == 0) && all(timing.destroyedInputSamples);
end
end

function freePlans(modules,plans)
for iPlan = 1:numel(plans)
    if plans(iPlan) ~= 0
        try
            feval(modules(iPlan),'free',plans(iPlan));
        catch
        end
    end
end
end

function comparison = failedComparison(configuration,nWarmups,nSamples,exception)
comparison = emptyComparison;
comparison.status = "failed";
comparison.configuration = configuration;
comparison.nWarmups = nWarmups;
comparison.nSamples = nSamples;
comparison.failure.identifier = string(exception.identifier);
comparison.failure.message = string(exception.message);
end

function workloads = selectWorkloadResults(comparisons,options)
workloads = repmat(emptyWorkload,size(options.sizes,1),1);
for iSize = 1:size(options.sizes,1)
    matches = arrayfun(@(comparison) comparison.status == "passed" && isequal(comparison.configuration.size,options.sizes(iSize,:)),comparisons);
    candidates = comparisons(matches);
    if isempty(candidates)
        error('FFTWBundledNative:NoValidComparison','No valid matched comparison remains for %s.',formatSize(options.sizes(iSize,:)));
    end
    validBundled = arrayfun(@(candidate) candidate.bundled.correctnessPassed && candidate.bundled.zeroCopyPassed,candidates);
    validNative = arrayfun(@(candidate) candidate.native.correctnessPassed && candidate.native.zeroCopyPassed,candidates);
    if ~any(validBundled) || ~any(validNative)
        error('FFTWBundledNative:NoValidEngineCandidate','No correct zero-copy candidate remains for %s.',formatSize(options.sizes(iSize,:)));
    end
    bundledIndex = rankedIndex(candidates,validBundled,"bundled");
    nativeIndex = rankedIndex(candidates,validNative,"native");
    selected = candidates(bundledIndex);
    workload = emptyWorkload;
    workload.size = options.sizes(iSize,:);
    workload.isGateWorkload = ismember(workload.size,options.gateSizes,'rows');
    workload.bestBundledConfiguration = selected.configuration;
    workload.bestBundled = selected.bundled;
    workload.matchedNative = selected.native;
    workload.bundledRelativeToMatchedNative = selected.bundledRelativeToNative;
    workload.bestNativeConfiguration = candidates(nativeIndex).configuration;
    workload.bestNative = candidates(nativeIndex).native;
    workload.matlab = selected.matlab;
    workload.storage = selected.storage;
    workload.thresholds.rawForwardPassed = workload.bestBundled.rawForwardSpeedRatio >= options.rawForwardSpeedThreshold;
    workload.thresholds.totalForwardPassed = workload.bestBundled.totalForwardSpeedRatio >= options.totalForwardSpeedThreshold;
    workload.thresholds.inversePassed = workload.bestBundled.inverseSpeedRatio >= options.inverseSpeedThreshold;
    workload.thresholds.correctnessPassed = workload.bestBundled.maximumRelativeError <= options.errorTolerance;
    workload.thresholds.storagePassed = workload.storage.halfSpectrumBytes < workload.storage.matlabFullSpectrumBytes;
    workload.thresholds.zeroCopyPassed = workload.bestBundled.zeroCopyPassed;
    workload.ready = workload.thresholds.rawForwardPassed && workload.thresholds.totalForwardPassed && workload.thresholds.inversePassed && workload.thresholds.correctnessPassed && workload.thresholds.storagePassed && workload.thresholds.zeroCopyPassed;
    workloads(iSize) = workload;
end
end

function index = rankedIndex(candidates,valid,engineField)
indices = find(valid);
total = arrayfun(@(candidate) candidate.(engineField).forward.totalMexMedianSeconds,candidates(indices));
raw = arrayfun(@(candidate) candidate.(engineField).forward.kernelMedianSeconds,candidates(indices));
[~,order] = sortrows([total(:) raw(:)],[1 2]);
index = indices(order(1));
end

function readiness = readinessRecord(workloads,options)
readiness = emptyReadiness;
gateMask = arrayfun(@(workload) ismember(workload.size,options.gateSizes,'rows'),workloads);
readiness.requiredGateWorkloads = size(options.gateSizes,1);
readiness.observedGateWorkloads = sum(gateMask);
readiness.rawForwardPassed = readiness.observedGateWorkloads == readiness.requiredGateWorkloads && all(arrayfun(@(workload) workload.thresholds.rawForwardPassed,workloads(gateMask)));
readiness.totalForwardPassed = readiness.observedGateWorkloads == readiness.requiredGateWorkloads && all(arrayfun(@(workload) workload.thresholds.totalForwardPassed,workloads(gateMask)));
readiness.inversePassed = readiness.observedGateWorkloads == readiness.requiredGateWorkloads && all(arrayfun(@(workload) workload.thresholds.inversePassed,workloads(gateMask)));
readiness.correctnessPassed = readiness.observedGateWorkloads == readiness.requiredGateWorkloads && all(arrayfun(@(workload) workload.thresholds.correctnessPassed,workloads(gateMask)));
readiness.storagePassed = readiness.observedGateWorkloads == readiness.requiredGateWorkloads && all(arrayfun(@(workload) workload.thresholds.storagePassed,workloads(gateMask)));
readiness.zeroCopyPassed = readiness.observedGateWorkloads == readiness.requiredGateWorkloads && all(arrayfun(@(workload) workload.thresholds.zeroCopyPassed,workloads(gateMask)));
readiness.ready = readiness.rawForwardPassed && readiness.totalForwardPassed && readiness.inversePassed && readiness.correctnessPassed && readiness.storagePassed && readiness.zeroCopyPassed;
if readiness.ready, readiness.status = "READY"; else, readiness.status = "NOT READY"; end
readiness.scope = "horizontal bundled-FFTW r2c/c2r only; issue #40 decides the complete support matrix";
end

function replays = replayIssue38(options,workloads)
if ~isfile(options.issue38ArtifactPath)
    error('FFTWBundledNative:Issue38ArtifactMissing','The issue #38 artifact was not found at %s.',options.issue38ArtifactPath);
end
artifact = jsondecode(fileread(options.issue38ArtifactPath));
candidates = normalizeStructArray(artifact.candidates);
oldWorkloads = normalizeStructArray(artifact.workloads);
replays = repmat(emptyHistoricalReplay,0,1);
for iSize = 1:size(options.sizes,1)
    sz = options.sizes(iSize,:);
    nSamples = options.nSamples;
    if iSize == size(options.sizes,1), nSamples = options.nSamplesLargest; end
    bundledMask = arrayfun(@(candidate) string(candidate.engine) == "bundled-fftw" && string(candidate.status) == "passed" && isequal(double(candidate.size(:).'),sz),candidates);
    bundledCandidates = candidates(bundledMask);
    if isempty(bundledCandidates), error('FFTWBundledNative:Issue38BundledCandidateMissing','No bundled issue #38 candidate matches %s.',formatSize(sz)); end
    screeningTimes = arrayfun(@(candidate) candidate.screening.forward.totalMexMedianSeconds,bundledCandidates);
    [~,bundledIndex] = min(screeningTimes);
    oldMatch = find(arrayfun(@(workload) isequal(double(workload.size(:).'),sz),oldWorkloads),1);
    if isempty(oldMatch), error('FFTWBundledNative:Issue38WinnerMissing','No issue #38 winner matches %s.',formatSize(sz)); end
    historical = {bundledCandidates(bundledIndex),oldWorkloads(oldMatch).winner};
    for iHistorical = 1:2
        candidate = historical{iHistorical};
        rng(options.seed+iSize-1,'twister');
        x = randn(sz);
        replay = replayHistoricalConfiguration(candidate,x,options.nWarmups,nSamples,options.plannerTimeLimitSeconds,options.errorTolerance);
        roles = ["best-bundled-screening","selected-native-winner"];
        replay.role = roles(iHistorical);
        replay.issue38ScreeningTotalMexSeconds = candidate.screening.forward.totalMexMedianSeconds;
        replay.issue38ScreeningRawSeconds = candidate.screening.forward.rawPipelineMedianSeconds;
        if iHistorical == 1
            replay.newMatchedTotalMexSeconds = workloads(iSize).bestBundled.forward.totalMexMedianSeconds;
        else
            replay.newMatchedTotalMexSeconds = workloads(iSize).bestNative.forward.totalMexMedianSeconds;
        end
        replay.diagnosticOnly = true;
        replays(end+1,1) = replay; %#ok<AGROW>
    end
end
end

function replay = replayHistoricalConfiguration(candidate,x,nWarmups,nSamples,timeLimit,errorTolerance)
engine = string(candidate.engine);
if engine == "bundled-fftw", module = "fftw_engine_benchmark_bundled"; else, module = "fftw_engine_benchmark_native"; end
configuration.size = double(candidate.size(:).');
configuration.transformOrder = double(candidate.transformOrder(:).');
configuration.layout = string(candidate.layout);
configuration.strategy = string(candidate.strategy);
configuration.planner = string(candidate.planner);
configuration.threads = double(candidate.threads);
configuration.alignmentMode = string(candidate.alignmentMode);
[referenceSpectrum,complexSize] = spectrumReference(x,configuration.transformOrder);
complexTemplate = complex(zeros(complexSize));
feval(module,'forgetWisdom');
[plan,~,scaleFactor,planningSeconds,limitReached,inputAlignment,outputAlignment] = feval(module,'create',configuration.size,configuration.transformOrder,char(configuration.strategy),configuration.threads,plannerFlags(configuration.planner),char(configuration.alignmentMode),timeLimit,x,complexTemplate);
planCleanup = onCleanup(@() feval(module,'free',plan));

operationIds = ["matlab-forward","matlab-inverse","engine-forward-allocating","engine-forward-preallocated","engine-inverse"];
total = nan(numel(operationIds),nSamples);
kernel = nan(numel(operationIds),nSamples);
pipeline = nan(numel(operationIds),nSamples);
matlabSpectrum = fftAlong(x,configuration.transformOrder);
matlabOutput = x;
allocated = complexTemplate;
preallocated = complexTemplate;
realOutput = zeros(configuration.size);
destroyed = referenceSpectrum;
for iRound = 1:(nWarmups+nSamples)
    firstOperation = mod(iRound-1,numel(operationIds))+1;
    order = [firstOperation:numel(operationIds) 1:firstOperation-1];
    for iOperation = order
        measuredKernel = NaN;
        measuredPipeline = NaN;
        switch operationIds(iOperation)
            case "matlab-forward"
                timer = tic; matlabSpectrum = fftAlong(x,configuration.transformOrder); elapsed = toc(timer);
            case "matlab-inverse"
                timer = tic; matlabOutput = ifftAlong(matlabSpectrum,configuration.transformOrder); elapsed = toc(timer);
            case "engine-forward-allocating"
                timer = tic; [allocated,measuredKernel,measuredPipeline] = feval(module,'r2c',plan,x); elapsed = toc(timer);
            case "engine-forward-preallocated"
                timer = tic; [preallocated,measuredKernel,measuredPipeline] = feval(module,'r2c',plan,x,preallocated); elapsed = toc(timer);
            case "engine-inverse"
                destroyed = complex(real(referenceSpectrum),imag(referenceSpectrum));
                timer = tic; [destroyed,realOutput,measuredKernel,measuredPipeline] = feval(module,'c2r',plan,destroyed,realOutput); elapsed = toc(timer);
        end
        if iRound > nWarmups
            iSample = iRound-nWarmups;
            total(iOperation,iSample) = elapsed;
            kernel(iOperation,iSample) = measuredKernel;
            pipeline(iOperation,iSample) = measuredPipeline;
        end
    end
end
allocatingMedian = median(total(3,:));
preallocatedMedian = median(total(4,:));
if allocatingMedian <= preallocatedMedian, forwardIndex = 3; allocationMode = "allocating"; else, forwardIndex = 4; allocationMode = "preallocated"; end
[~,allocatingError] = numericalErrors(allocated,referenceSpectrum);
[~,preallocatedError] = numericalErrors(preallocated,referenceSpectrum);
[~,inverseError] = numericalErrors(scaleFactor*realOutput,x);
replay = emptyHistoricalReplay;
replay.status = "passed";
replay.engine = engine;
replay.configuration = configuration;
replay.nWarmups = nWarmups;
replay.nSamples = nSamples;
replay.plan.planningSeconds = planningSeconds;
replay.plan.planningLimitReached = logical(limitReached);
replay.plan.inputAlignmentClass = inputAlignment;
replay.plan.outputAlignmentClass = outputAlignment;
replay.plan.wisdomClearedBeforePlanning = true;
replay.matlabForwardMedianSeconds = median(total(1,:));
replay.matlabInverseMedianSeconds = median(total(2,:));
replay.forwardAllocationMode = allocationMode;
replay.forwardTotalMexSamplesSeconds = total(forwardIndex,:);
replay.forwardTotalMexMedianSeconds = median(total(forwardIndex,:));
replay.forwardKernelSamplesSeconds = kernel(forwardIndex,:);
replay.forwardKernelMedianSeconds = median(kernel(forwardIndex,:));
replay.forwardRawPipelineSamplesSeconds = pipeline(forwardIndex,:);
replay.forwardRawPipelineMedianSeconds = median(pipeline(forwardIndex,:));
replay.inverseTotalMexSamplesSeconds = total(5,:);
replay.inverseTotalMexMedianSeconds = median(total(5,:));
replay.maximumRelativeError = max([allocatingError preallocatedError inverseError numericalError(matlabOutput,x)]);
replay.correctnessPassed = replay.maximumRelativeError <= errorTolerance;
replay.destructiveInputChanged = ~isequaln(destroyed,referenceSpectrum);
clear planCleanup
end

function values = normalizeStructArray(values)
if iscell(values), values = [values{:}]; end
values = values(:);
end

function flags = plannerFlags(planner)
switch planner
    case "estimate", flags = 64;
    case "measure", flags = 0;
    case "patient", flags = 32;
    case "exhaustive", flags = 8;
    otherwise, error('FFTWBundledNative:UnknownPlanner','Unknown planner %s.',planner);
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

function relativeError = numericalError(actual,reference)
[~,relativeError] = numericalErrors(actual,reference);
end

function [absoluteError,relativeError] = numericalErrors(actual,reference)
absoluteError = max(abs(actual-reference),[],'all');
scale = max(abs(reference),[],'all');
if scale == 0, scale = 1; end
relativeError = absoluteError/scale;
end

function storage = storageRecord(realSize,complexSize)
storage.realInputBytes = prod(realSize)*8;
storage.matlabFullSpectrumBytes = prod(realSize)*16;
storage.halfSpectrumBytes = prod(complexSize)*16;
storage.realOutputBytes = prod(realSize)*8;
storage.halfToFullSpectrumRatio = storage.halfSpectrumBytes/storage.matlabFullSpectrumBytes;
end

function value = percentDifference(value,reference)
value = 100*(value/reference-1);
end

function environment = collectEnvironment()
environment.matlabVersion = string(version);
environment.matlabRelease = string(version('-release'));
environment.architecture = string(computer('arch'));
environment.mexExtension = string(mexext);
environment.operatingSystem = string(system_dependent('getos'));
environment.processor = "unknown";
environment.machineModel = "unknown";
environment.machineName = "unknown";
environment.physicalMemory = "unknown";
[status,hardwareText] = system('system_profiler SPHardwareDataType -json');
if status == 0
    data = jsondecode(hardwareText);
    hardware = data.SPHardwareDataType(1);
    environment.processor = optionalStructString(hardware,'chip_type',environment.processor);
    environment.machineModel = optionalStructString(hardware,'machine_model',environment.machineModel);
    environment.machineName = optionalStructString(hardware,'machine_name',environment.machineName);
    environment.physicalMemory = optionalStructString(hardware,'physical_memory',environment.physicalMemory);
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
if isfield(inputStruct,fieldName), value = string(inputStruct.(fieldName)); else, value = defaultValue; end
end

function sources = collectSourceProvenance(sourceDirectory)
names = ["runFFTWBundledNativeComparison.m","buildFFTWBundledNativeComparisonMex.m","buildFFTWMexOwnershipBenchmark.m","fftw_ownership_benchmark.cpp","runFFTWEngineLayoutBenchmark.m","fftw_engine_benchmark.cpp"];
sources.files = repmat(struct('path',"",'sha256',""),numel(names),1);
for iFile = 1:numel(names)
    sources.files(iFile).path = names(iFile);
    sources.files(iFile).sha256 = fileSHA256(fullfile(sourceDirectory,names(iFile)));
end
repositoryRoot = fileparts(fileparts(fileparts(sourceDirectory)));
[status,commit] = system(sprintf('git -C "%s" rev-parse HEAD',repositoryRoot));
if status == 0, sources.repositoryCommit = string(strtrim(commit)); else, sources.repositoryCommit = "unknown"; end
[status,dirty] = system(sprintf('git -C "%s" status --porcelain --untracked-files=all',repositoryRoot));
sources.repositoryDirty = status ~= 0 || strlength(strtrim(string(dirty))) > 0;
end

function hash = fileSHA256(path)
digest = java.security.MessageDigest.getInstance('SHA-256');
digest.update(unicode2native(fileread(path),'UTF-8'));
hashBytes = typecast(digest.digest(),'uint8');
hash = string(lower(reshape(dec2hex(hashBytes,2).',1,[])));
end

function writeResultArtifacts(result,runDirectory)
serializable = replaceMissingStrings(result);
writeTextFile(fullfile(runDirectory,"bundled-native-comparison.json"),string(jsonencode(serializable)));
writeTextFile(fullfile(runDirectory,"summary.md"),createMarkdownSummary(result));
end

function value = replaceMissingStrings(value)
if isstring(value)
    value(ismissing(value)) = "";
elseif isstruct(value)
    fields = fieldnames(value);
    for iValue = 1:numel(value)
        for iField = 1:numel(fields)
            value(iValue).(fields{iField}) = replaceMissingStrings(value(iValue).(fields{iField}));
        end
    end
elseif iscell(value)
    for iValue = 1:numel(value)
        value{iValue} = replaceMissingStrings(value{iValue});
    end
end
end

function summary = createMarkdownSummary(result)
lines = strings(0,1);
lines(end+1) = "# Bundled versus native FFTW comparison";
lines(end+1) = "";
lines(end+1) = "Benchmark status: **" + upper(result.status) + "**";
if result.status == "passed", lines(end+1) = "Horizontal bundled-FFTW readiness: **" + result.readiness.status + "**"; end
lines(end+1) = "";
lines(end+1) = "Native FFTW is a controlled reference only and is not a runtime requirement.";
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
if isfield(result.configuration,'observedThreads'), lines(end+1) = markdownRow("Threads",string(result.configuration.observedThreads)); end
lines(end+1) = "";
lines(end+1) = "## Engines";
lines(end+1) = "";
lines(end+1) = "| Engine | Version | Resolved library | Matched alignment | Mismatch rejected | Unaligned |";
lines(end+1) = "|---|---|---|:---:|:---:|:---:|";
for engine = result.engines'
    lines(end+1) = sprintf('| %s | %s | %s | %s | %s | %s |',engine.id,escapeMarkdown(engine.version),escapeMarkdown(engine.library),yesNo(engine.alignmentMatchedAccepted),yesNo(engine.alignmentMismatchRejected),yesNo(engine.unalignedAccepted)); %#ok<AGROW>
end
if ~isempty(result.historicalReplay)
    lines(end+1) = "";
    lines(end+1) = "## Issue #38 discrepancy";
    lines(end+1) = "";
    lines(end+1) = "Issue #38 selected one candidate from a large matrix using three screening samples, recreated that plan without an explicit wisdom reset, changed from a three-operation to a five-operation schedule, and did not persist the bundled finalist timings. The rows below replay the reconstructable finalists with isolated wisdom and full sampling; staged results are diagnostic only.";
    lines(end+1) = "";
    lines(end+1) = "| Size | Role | Engine | Strategy | Planner | Alignment | #38 screen (ms) | Full replay (ms) | New matched result (ms) |";
    lines(end+1) = "|---|---|---|---|---|---|---:|---:|---:|";
    for replay = result.historicalReplay'
        lines(end+1) = sprintf('| %s | %s | %s | %s | %s | %s | %.3f | %.3f | %.3f |',formatSize(replay.configuration.size),replay.role,replay.engine,replay.configuration.strategy,replay.configuration.planner,replay.configuration.alignmentMode,1e3*replay.issue38ScreeningTotalMexSeconds,1e3*replay.forwardTotalMexMedianSeconds,1e3*replay.newMatchedTotalMexSeconds); %#ok<AGROW>
    end
end
if ~isempty(result.comparisons)
    lines(end+1) = "";
    lines(end+1) = "## Matched configuration comparison";
    lines(end+1) = "";
    lines(end+1) = "`Bundled vs native total` is the percentage by which bundled time differs from native time; positive means bundled FFTW was slower.";
    lines(end+1) = "";
    lines(end+1) = "| Size | Planner | Alignment | Status | MATLAB forward (ms) | Bundled raw (ms) | Bundled total (ms) | Native total (ms) | Bundled vs native total |";
    lines(end+1) = "|---|---|---|---|---:|---:|---:|---:|---:|";
    for comparison = result.comparisons'
        if comparison.status == "passed"
            lines(end+1) = sprintf('| %s | %s | %s | %s | %.3f | %.3f | %.3f | %.3f | %+.1f%% |',formatSize(comparison.configuration.size),comparison.configuration.planner,comparison.configuration.alignmentMode,comparison.status,1e3*comparison.matlab.forward.medianSeconds,1e3*comparison.bundled.forward.kernelMedianSeconds,1e3*comparison.bundled.forward.totalMexMedianSeconds,1e3*comparison.native.forward.totalMexMedianSeconds,comparison.bundledRelativeToNative.forwardTotalPercent); %#ok<AGROW>
        else
            lines(end+1) = sprintf('| %s | %s | %s | failed | — | — | — | — | — |',formatSize(comparison.configuration.size),comparison.configuration.planner,comparison.configuration.alignmentMode); %#ok<AGROW>
        end
    end
end
if ~isempty(result.workloads)
    lines(end+1) = "";
    lines(end+1) = "## Readiness by workload";
    lines(end+1) = "";
    lines(end+1) = "Speedups above 1 are faster than MATLAB. Positive bundled/native percentages mean bundled FFTW was slower.";
    lines(end+1) = "";
    lines(end+1) = "| Size | Best bundled configuration | MATLAB forward (ms) | Bundled raw (ms) | Raw speedup | Bundled total (ms) | Total speedup | Bundled time relative to native | Destructive c2r speedup | Relative error | Half/full storage | Result |";
    lines(end+1) = "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|";
    for workload = result.workloads'
        configuration = workload.bestBundledConfiguration;
        lines(end+1) = sprintf('| %s | %s / %s / %d threads | %.3f | %.3f | %.3fx | %.3f | %.3fx | %+.1f%% | %.3fx | %.3g | %.3f | %s |',formatSize(workload.size),configuration.planner,configuration.alignmentMode,configuration.threads,1e3*workload.matlab.forward.medianSeconds,1e3*workload.bestBundled.forward.kernelMedianSeconds,workload.bestBundled.rawForwardSpeedRatio,1e3*workload.bestBundled.forward.totalMexMedianSeconds,workload.bestBundled.totalForwardSpeedRatio,workload.bundledRelativeToMatchedNative.forwardTotalPercent,workload.bestBundled.inverseSpeedRatio,workload.bestBundled.maximumRelativeError,workload.storage.halfToFullSpectrumRatio,readyText(workload.ready)); %#ok<AGROW>
    end
    lines(end+1) = "";
    lines(end+1) = "## Criterion status";
    lines(end+1) = "";
    lines(end+1) = "| Size | Raw r2c >= 1.25x | Total MEX r2c >= 1.10x | Destructive c2r >= 0.95x | Error <= 1e-12 | Half spectrum | Zero copy |";
    lines(end+1) = "|---|:---:|:---:|:---:|:---:|:---:|:---:|";
    for workload = result.workloads'
        threshold = workload.thresholds;
        lines(end+1) = sprintf('| %s | %s | %s | %s | %s | %s | %s |',formatSize(workload.size),yesNo(threshold.rawForwardPassed),yesNo(threshold.totalForwardPassed),yesNo(threshold.inversePassed),yesNo(threshold.correctnessPassed),yesNo(threshold.storagePassed),yesNo(threshold.zeroCopyPassed)); %#ok<AGROW>
    end
end
lines(end+1) = "";
lines(end+1) = "## Timing boundaries";
lines(end+1) = "";
lines(end+1) = "- Raw r2c time is the FFTW execute call only; no packing or layout conversion is required.";
lines(end+1) = "- Internal pipeline time includes MATLAB-buffer allocation, FFT execution, and buffer wrapping.";
lines(end+1) = "- Complete MEX time includes the entire MATLAB-to-MEX call. Input generation, planning, destructive-input refresh, metrics retrieval, and inverse normalization are excluded.";
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

function writeTextFile(path,text)
[fileId,message] = fopen(path,'w');
if fileId < 0, error('FFTWBundledNative:ArtifactOpenFailed','Unable to open %s: %s',path,message); end
cleanup = onCleanup(@() fclose(fileId));
fprintf(fileId,'%s\n',text);
clear cleanup
end

function restoreMatlabState(planner,wisdom,threads,randomState)
try, fftw('dwisdom',[]); fftw('dwisdom',wisdom); catch, end
try, fftw('planner',char(planner)); catch, end
try, maxNumCompThreads(threads); catch, end
try, rng(randomState); catch, end
end

function failure = failureRecord(exception,stage)
failure.stage = stage;
failure.identifier = string(exception.identifier);
failure.message = string(exception.message);
failure.stack = arrayfun(@(entry) struct('file',string(entry.file),'name',string(entry.name),'line',entry.line),exception.stack);
end

function value = utcTimestamp()
value = string(datetime('now',TimeZone='UTC',Format="yyyy-MM-dd'T'HH:mm:ss.SSS'Z'"));
end

function value = layoutName(compressedDimension)
names = ["half-x","half-y","half-z"];
value = names(compressedDimension);
end

function value = formatSize(sz)
value = strtrim(sprintf('%d x %d x %d',sz));
end

function row = markdownRow(field,value)
row = "| " + escapeMarkdown(string(field)) + " | " + escapeMarkdown(string(value)) + " |";
end

function value = escapeMarkdown(value)
value = replace(string(value),"|","\|");
value = replace(value,newline,"<br>");
end

function value = yesNo(condition)
if condition, value = "yes"; else, value = "no"; end
end

function value = readyText(condition)
if condition, value = "READY"; else, value = "NOT READY"; end
end

function engine = emptyEngine()
engine = struct('id',"",'module',"",'version',"",'library',"",'alignmentMatchedAccepted',false,'alignmentMismatchRejected',false,'unalignedAccepted',false);
end

function plan = emptyPlan()
plan = struct('planningSeconds',NaN,'planningLimitReached',false,'inputAlignmentClass',NaN,'outputAlignmentClass',NaN,'wisdomClearedBeforePlanning',false);
end

function timing = emptyEngineTiming()
timing = struct('totalMexSamplesSeconds',[],'totalMexMedianSeconds',NaN,'allocationSamplesSeconds',[],'allocationMedianSeconds',NaN,'bufferWrapSamplesSeconds',[],'bufferWrapMedianSeconds',NaN,'explicitMemcpySamplesSeconds',[],'explicitMemcpyMedianSeconds',NaN,'mutableDetachSamplesSeconds',[],'mutableDetachMedianSeconds',NaN,'kernelSamplesSeconds',[],'kernelMedianSeconds',NaN,'internalPipelineSamplesSeconds',[],'internalPipelineMedianSeconds',NaN,'boundaryResidualSamplesSeconds',[],'boundaryResidualMedianSeconds',NaN,'explicitCopyCountSamples',[],'explicitCopiedBytesSamples',[],'detectedCopyCountSamples',[],'detectedCopiedBytesSamples',[],'inputAlignmentClasses',[],'outputAlignmentClasses',[],'destroyedInputSamples',[],'inputBeforeTokens',uint64([]),'inputMutableTokens',uint64([]),'outputMutableTokens',uint64([]),'wrappedTokens',uint64([]),'returnedTokens',uint64([]),'pointerPreserved',false,'zeroCopyPassed',false,'maximumAbsoluteError',NaN,'maximumRelativeError',NaN,'correctnessPassed',false);
end

function measurement = emptyEngineMeasurement()
measurement = struct('plan',emptyPlan,'forward',emptyEngineTiming,'inverse',emptyEngineTiming,'maximumRelativeError',NaN,'correctnessPassed',false,'zeroCopyPassed',false,'destructiveInputChanged',false,'rawForwardSpeedRatio',NaN,'totalForwardSpeedRatio',NaN,'inverseSpeedRatio',NaN);
end

function comparison = emptyComparison()
comparison = struct('status',"",'configuration',struct,'nWarmups',0,'nSamples',0,'matlab',struct('forward',struct('samplesSeconds',[],'medianSeconds',NaN),'inverse',struct('samplesSeconds',[],'medianSeconds',NaN),'maximumRelativeError',NaN),'bundled',emptyEngineMeasurement,'native',emptyEngineMeasurement,'bundledRelativeToNative',struct('forwardKernelPercent',NaN,'forwardTotalPercent',NaN,'inverseTotalPercent',NaN),'storage',struct,'failure',struct('identifier',"",'message',""));
end

function workload = emptyWorkload()
thresholds = struct('rawForwardPassed',false,'totalForwardPassed',false,'inversePassed',false,'correctnessPassed',false,'storagePassed',false,'zeroCopyPassed',false);
workload = struct('size',[],'isGateWorkload',false,'bestBundledConfiguration',struct,'bestBundled',emptyEngineMeasurement,'matchedNative',emptyEngineMeasurement,'bundledRelativeToMatchedNative',struct,'bestNativeConfiguration',struct,'bestNative',emptyEngineMeasurement,'matlab',struct,'storage',struct,'thresholds',thresholds,'ready',false);
end

function replay = emptyHistoricalReplay()
replay = struct('status',"",'role',"",'engine',"",'configuration',struct,'nWarmups',0,'nSamples',0,'plan',emptyPlan,'matlabForwardMedianSeconds',NaN,'matlabInverseMedianSeconds',NaN,'forwardAllocationMode',"",'forwardTotalMexSamplesSeconds',[],'forwardTotalMexMedianSeconds',NaN,'forwardKernelSamplesSeconds',[],'forwardKernelMedianSeconds',NaN,'forwardRawPipelineSamplesSeconds',[],'forwardRawPipelineMedianSeconds',NaN,'inverseTotalMexSamplesSeconds',[],'inverseTotalMexMedianSeconds',NaN,'maximumRelativeError',NaN,'correctnessPassed',false,'destructiveInputChanged',false,'issue38ScreeningTotalMexSeconds',NaN,'issue38ScreeningRawSeconds',NaN,'newMatchedTotalMexSeconds',NaN,'diagnosticOnly',true);
end

function readiness = emptyReadiness()
readiness = struct('status',"NOT READY",'ready',false,'requiredGateWorkloads',0,'observedGateWorkloads',0,'rawForwardPassed',false,'totalForwardPassed',false,'inversePassed',false,'correctnessPassed',false,'storagePassed',false,'zeroCopyPassed',false,'scope',"");
end
