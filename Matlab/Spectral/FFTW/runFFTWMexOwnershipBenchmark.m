function result = runFFTWMexOwnershipBenchmark(options)
% Characterize FFTW MEX ownership, allocation, and copy costs.
%
% This experimental issue #39 benchmark imports the issue #38 winning
% configurations and measures MATLAB/MEX ownership behavior. It does not
% make the milestone's GO or NO-GO decision.
arguments
    options.outputDirectory (1,1) string = fullfile(fileparts(mfilename('fullpath')),"results","issue39")
    options.issue38ArtifactPath (1,1) string = fullfile(fileparts(mfilename('fullpath')),"results","issue38","20260808T045405991Z-maca64-r2026a","engine-layout-benchmark.json")
    options.sizes (:,3) double = [128 128 128; 256 128 40; 256 256 64; 1024 1024 30]
    options.gateSizes (:,3) double = [256 256 64; 1024 1024 30]
    options.nWarmups (1,1) double = 2
    options.nSamples (1,1) double = 7
    options.nSamplesLargest (1,1) double = 3
    options.nPlanMemorySamples (1,1) double = 3
    options.seed (1,1) double = 37
    options.errorTolerance (1,1) double = 1e-12
    options.shouldUseIssue38Winners (1,1) logical = true
    options.smokeEngine (1,1) string = "bundled-fftw"
    options.smokeTransformOrder (1,2) double = [2 1]
    options.smokePlanner (1,1) string = "estimate"
    options.smokeThreads (1,1) double = 1
    options.smokeAlignmentMode (1,1) string = "unaligned"
    options.shouldBuild (1,1) logical = true
    options.shouldBuildNative (1,1) logical = true
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
        result.build = buildFFTWMexOwnershipBenchmark(shouldBuildNative=options.shouldBuildNative);
    else
        requireModule("fftw_ownership_benchmark_bundled");
        result.build = existingBuildRecord();
    end
    result.engines = collectEngineRecords(result.build);
    configurations = ownershipConfigurations(options,result.build);
    result.configuration.workloadConfigurations = configurations;

    activeStage = "ownership measurements";
    for iSize = 1:size(options.sizes,1)
        nSamples = options.nSamples;
        if iSize == size(options.sizes,1), nSamples = options.nSamplesLargest; end
        activeStage = "ownership measurement " + formatSize(options.sizes(iSize,:));
        result.workloads(end+1,1) = benchmarkWorkload(configurations(iSize),options.sizes(iSize,:),options.seed+iSize-1,options.nWarmups,nSamples,options.nPlanMemorySamples,options.errorTolerance);
    end
    result.recommendation = chooseRecommendation(result.workloads,options.gateSizes,options.errorTolerance);
    result.status = "passed";
    result.completedAtUTC = utcTimestamp();
    if options.shouldWriteArtifacts, writeResultArtifacts(result,runDirectory); end
catch exception
    result.status = "failed";
    result.completedAtUTC = utcTimestamp();
    result.failure = failureRecord(exception,activeStage);
    if options.shouldWriteArtifacts
        try
            writeResultArtifacts(result,runDirectory);
        catch artifactException
            warning('FFTWMexOwnership:FailureArtifactWriteFailed','Unable to write issue #39 failure artifacts: %s',artifactException.message);
        end
    end
    rethrow(exception);
end
clear stateCleanup
end

function validateOptions(options)
validateattributes(options.sizes,{'double'},{'2d','ncols',3,'integer','positive'},mfilename,'sizes');
validateattributes(options.gateSizes,{'double'},{'2d','ncols',3,'integer','positive'},mfilename,'gateSizes');
validateattributes(options.nWarmups,{'double'},{'scalar','integer','nonnegative'},mfilename,'nWarmups');
validateattributes(options.nSamples,{'double'},{'scalar','integer','positive'},mfilename,'nSamples');
validateattributes(options.nSamplesLargest,{'double'},{'scalar','integer','positive'},mfilename,'nSamplesLargest');
validateattributes(options.nPlanMemorySamples,{'double'},{'scalar','integer','positive'},mfilename,'nPlanMemorySamples');
validateattributes(options.seed,{'double'},{'scalar','integer','nonnegative'},mfilename,'seed');
validateattributes(options.errorTolerance,{'double'},{'scalar','real','finite','nonnegative'},mfilename,'errorTolerance');
validateattributes(options.smokeThreads,{'double'},{'scalar','integer','positive'},mfilename,'smokeThreads');
if options.smokeTransformOrder(1) == options.smokeTransformOrder(2) || any(options.smokeTransformOrder > 3)
    error('FFTWMexOwnership:InvalidTransformOrder','The transform order must contain two distinct dimensions between 1 and 3.');
end
if ~ismember(options.smokeEngine,["bundled-fftw","native-fftw"])
    error('FFTWMexOwnership:UnknownEngine','The smoke engine must be bundled-fftw or native-fftw.');
end
if ~ismember(options.smokePlanner,["estimate","measure","patient","exhaustive"])
    error('FFTWMexOwnership:UnknownPlanner','Unknown FFTW planner.');
end
if ~ismember(options.smokeAlignmentMode,["matched","unaligned"])
    error('FFTWMexOwnership:UnknownAlignmentMode','Alignment mode must be matched or unaligned.');
end
if ~startsWith(string(version('-release')),"2026a",IgnoreCase=true)
    error('FFTWMexOwnership:UnsupportedMATLABRelease','Issue #39 targets MATLAB R2026a only.');
end
end

function runDirectory = prepareRunDirectory(outputDirectory,runId,shouldWrite)
runDirectory = "";
if ~shouldWrite, return; end
runDirectory = fullfile(outputDirectory,runId);
if isfolder(runDirectory) || isfile(runDirectory)
    error('FFTWMexOwnership:OutputExists','Benchmark output already exists at %s.',runDirectory);
end
[didCreate,message] = mkdir(runDirectory);
if ~didCreate, error('FFTWMexOwnership:OutputCreationFailed','Unable to create %s: %s',runDirectory,message); end
end

function result = initializeResult(runId,runDirectory,sourceDirectory,options)
result.schemaVersion = "1.0.0";
result.status = "running";
result.runId = runId;
result.generatedAtUTC = utcTimestamp();
result.completedAtUTC = "";
result.environment = collectEnvironment();
result.configuration.sizes = options.sizes;
result.configuration.gateSizes = options.gateSizes;
result.configuration.issue38ArtifactPath = options.issue38ArtifactPath;
result.configuration.nWarmups = options.nWarmups;
result.configuration.nSamples = options.nSamples;
result.configuration.nSamplesLargest = options.nSamplesLargest;
result.configuration.nPlanMemorySamples = options.nPlanMemorySamples;
result.configuration.seeds = options.seed + (0:size(options.sizes,1)-1);
result.configuration.errorTolerance = options.errorTolerance;
result.configuration.matlabCompatibility = "R2026a only";
result.configuration.operationOrder = "round-robin with rotating first operation";
result.configuration.timingBoundary = "complete MEX call; metrics retrieval excluded";
result.configuration.recommendationRule = "valid zero-copy gate candidates ranked by geometric-mean complete MEX time; if common FFT-kernel noise leaves no model within 5 percent on both gates, use common-kernel plus measured ownership overhead; within 3 percent prefer direct caller reuse, then other MATLAB ownership, then FFTW custom ownership";
result.sources = collectSourceProvenance(sourceDirectory);
result.build = struct;
result.engines = repmat(emptyEngine(),0,1);
result.workloads = repmat(emptyWorkload(),0,1);
result.recommendation = emptyRecommendation();
result.failure = [];
result.artifacts.directory = runDirectory;
result.artifacts.json = "ownership-benchmark.json";
result.artifacts.markdown = "summary.md";
end

function build = existingBuildRecord()
build.bundled.module = "fftw_ownership_benchmark_bundled";
build.bundled.status = "existing";
build.bundled.library = string(fullfile(matlabroot,"bin",computer('arch'),"libmwfftw3.3.dylib"));
build.native.module = "fftw_ownership_benchmark_native";
if exist('fftw_ownership_benchmark_native','file') == 3
    build.native.status = "existing";
    build.native.reason = "";
else
    build.native.status = "unavailable";
    build.native.reason = "The native ownership module is not built.";
end
end

function requireModule(module)
if exist(module,'file') ~= 3, error('FFTWMexOwnership:MexMissing','%s is unavailable and shouldBuild is false.',module); end
end

function engines = collectEngineRecords(build)
engines = repmat(emptyEngine(),0,1);
modules = ["fftw_ownership_benchmark_bundled","fftw_ownership_benchmark_native"];
ids = ["bundled-fftw","native-fftw"];
for iEngine = 1:numel(modules)
    if exist(modules(iEngine),'file') == 3
        [engineVersion,library] = feval(modules(iEngine),'info');
        engine = emptyEngine();
        engine.id = ids(iEngine);
        engine.module = modules(iEngine);
        engine.status = "available";
        engine.version = string(engineVersion);
        engine.library = string(library);
        [engine.alignmentMatchedAccepted,engine.alignmentMismatchRejected,engine.unalignedAccepted] = feval(modules(iEngine),'alignmentSelfTest');
        engines(end+1,1) = engine; %#ok<AGROW>
    end
end
if isempty(engines), error('FFTWMexOwnership:NoEngine','No ownership benchmark FFTW engine is available.'); end
if isfield(build,'native') && isfield(build.native,'status') && string(build.native.status) == "built"
    nativeIndex = find([engines.id] == "native-fftw",1);
    if ~isempty(nativeIndex) && startsWith(engines(nativeIndex).library,string(matlabroot))
        error('FFTWMexOwnership:NativeResolutionFailure','The native module resolved FFTW symbols to MATLAB''s bundled library.');
    end
end
end

function configurations = ownershipConfigurations(options,build)
configurations = repmat(emptyConfiguration(),size(options.sizes,1),1);
if options.shouldUseIssue38Winners
    if ~isfile(options.issue38ArtifactPath), error('FFTWMexOwnership:Issue38ArtifactMissing','Issue #38 artifact not found at %s.',options.issue38ArtifactPath); end
    issue38 = jsondecode(fileread(options.issue38ArtifactPath));
    for iSize = 1:size(options.sizes,1)
        match = find(arrayfun(@(workload) isequal(double(workload.size(:).'),options.sizes(iSize,:)),issue38.workloads),1);
        if isempty(match), error('FFTWMexOwnership:Issue38WinnerMissing','No issue #38 winner matches %s.',formatSize(options.sizes(iSize,:))); end
        winner = issue38.workloads(match).winner;
        configurations(iSize) = configurationRecord(options.sizes(iSize,:),string(winner.engine),double(winner.transformOrder(:).'),string(winner.planner),double(winner.threads),string(winner.alignmentMode));
    end
else
    for iSize = 1:size(options.sizes,1)
        configurations(iSize) = configurationRecord(options.sizes(iSize,:),options.smokeEngine,options.smokeTransformOrder,options.smokePlanner,min(options.smokeThreads,maxNumCompThreads),options.smokeAlignmentMode);
    end
end
for iConfiguration = 1:numel(configurations)
    if configurations(iConfiguration).engine == "bundled-fftw"
        configurations(iConfiguration).module = "fftw_ownership_benchmark_bundled";
    else
        configurations(iConfiguration).module = "fftw_ownership_benchmark_native";
        if ~isfield(build,'native') || ~ismember(string(build.native.status),["built","existing"])
            error('FFTWMexOwnership:NativeModuleMissing','The issue #38 winner requires native FFTW, but the ownership module is unavailable.');
        end
    end
    requireModule(configurations(iConfiguration).module);
end
end

function configuration = configurationRecord(sz,engine,transformOrder,planner,threads,alignmentMode)
configuration = emptyConfiguration();
configuration.size = sz;
configuration.engine = engine;
configuration.transformOrder = transformOrder;
configuration.layout = "half-" + string(char('x'+transformOrder(2)-1));
configuration.strategy = "guru-rank2";
configuration.planner = planner;
configuration.threads = threads;
configuration.alignmentMode = alignmentMode;
end

function workload = benchmarkWorkload(configuration,sz,seed,nWarmups,nSamples,nPlanMemorySamples,errorTolerance)
module = configuration.module;
rng(seed,'twister');
x = randn(sz);
fullReference = fftAlong(x,configuration.transformOrder);
complexSize = sz;
compressedDimension = configuration.transformOrder(2);
complexSize(compressedDimension) = floor(sz(compressedDimension)/2)+1;
subscripts = repmat({':'},1,numel(sz));
subscripts{compressedDimension} = 1:complexSize(compressedDimension);
referenceSpectrum = fullReference(subscripts{:});
realTemplate = zeros(sz);
complexTemplate = complex(zeros(complexSize));
planner = plannerFlags(configuration.planner);

planMemory = probePlanMemory(module,sz,configuration.transformOrder,configuration.threads,planner,configuration.alignmentMode,realTemplate,complexTemplate,nPlanMemorySamples,numel(referenceSpectrum)*16);
[plan,reportedSize,scaleFactor,planningSeconds,inputAlignment,outputAlignment] = feval(module,'create',sz,configuration.transformOrder,configuration.threads,planner,char(configuration.alignmentMode),realTemplate,complexTemplate);
planCleanup = onCleanup(@() feval(module,'free',plan));
if ~isequal(double(reportedSize),double(complexSize)), error('FFTWMexOwnership:ReportedSizeMismatch','The ownership gateway reported an unexpected spectrum size.'); end
feval(module,'resetLifetime');

operationIds = ["forward-factory-array","forward-caller-direct","forward-caller-complex","forward-matlab-buffer","forward-fftw-buffer","inverse-preserving-allocating","inverse-preserving-preallocated","inverse-destructive-unique","inverse-destructive-aliased","inverse-destructive-complex","chain-fftw-forward-destructive-inverse"];
nOperations = numel(operationIds);
totalSamples = nan(nOperations,nSamples);
metricSamples = nan(nOperations,17,nSamples);
pointerSamples = zeros(nOperations,5,nSamples,'uint64');
returnedTokens = zeros(nOperations,nSamples,'uint64');
wrapperOriginalTokens = zeros(nOperations,nSamples,'uint64');

forwardFactory = complexTemplate;
forwardCaller = complexTemplate;
forwardComplex = complexTemplate;
forwardMatlabBuffer = complexTemplate;
forwardFFTWBuffer = complexTemplate;
inverseAllocating = realTemplate;
inversePreallocated = realTemplate;
inverseUniqueOutput = realTemplate;
inverseAliasedOutput = realTemplate;
inverseComplexOutput = realTemplate;
chainOutput = realTemplate;
preservingInput = complex(real(referenceSpectrum),imag(referenceSpectrum));
destroyedUnique = preservingInput;
aliasedInput = preservingInput;
chainSpectrum = complexTemplate;

for iRound = 1:(nWarmups+nSamples)
    firstOperation = mod(iRound-1,nOperations)+1;
    operationOrder = [firstOperation:nOperations 1:firstOperation-1];
    for iOperation = operationOrder
        wrapperOriginal = uint64(0);
        switch operationIds(iOperation)
            case "forward-factory-array"
                timer = tic;
                forwardFactory = feval(module,'forward',plan,x,'factory-array');
                elapsed = toc(timer);
                [metrics,pointers] = feval(module,'metrics',plan);
                returned = feval(module,'pointer',forwardFactory);
            case "forward-caller-direct"
                timer = tic;
                forwardCaller = feval(module,'forward',plan,x,'caller-direct',forwardCaller);
                elapsed = toc(timer);
                [metrics,pointers] = feval(module,'metrics',plan);
                returned = feval(module,'pointer',forwardCaller);
            case "forward-caller-complex"
                wrapperOriginal = feval(module,'pointer',forwardComplex);
                timer = tic;
                forwardComplex = feval(module,'forward',plan,x,'caller-direct',complex(forwardComplex));
                elapsed = toc(timer);
                [metrics,pointers] = feval(module,'metrics',plan);
                returned = feval(module,'pointer',forwardComplex);
            case "forward-matlab-buffer"
                timer = tic;
                forwardMatlabBuffer = feval(module,'forward',plan,x,'matlab-buffer');
                elapsed = toc(timer);
                [metrics,pointers] = feval(module,'metrics',plan);
                returned = feval(module,'pointer',forwardMatlabBuffer);
            case "forward-fftw-buffer"
                timer = tic;
                forwardFFTWBuffer = feval(module,'forward',plan,x,'fftw-buffer');
                elapsed = toc(timer);
                [metrics,pointers] = feval(module,'metrics',plan);
                returned = feval(module,'pointer',forwardFFTWBuffer);
            case "inverse-preserving-allocating"
                timer = tic;
                inverseAllocating = feval(module,'inversePreserving',plan,preservingInput,'allocating');
                elapsed = toc(timer);
                [metrics,pointers] = feval(module,'metrics',plan);
                returned = feval(module,'pointer',inverseAllocating);
            case "inverse-preserving-preallocated"
                timer = tic;
                inversePreallocated = feval(module,'inversePreserving',plan,preservingInput,'preallocated',inversePreallocated);
                elapsed = toc(timer);
                [metrics,pointers] = feval(module,'metrics',plan);
                returned = feval(module,'pointer',inversePreallocated);
            case "inverse-destructive-unique"
                destroyedUnique = complex(real(referenceSpectrum),imag(referenceSpectrum));
                timer = tic;
                [destroyedUnique,inverseUniqueOutput] = feval(module,'inverseDestructive',plan,destroyedUnique,inverseUniqueOutput);
                elapsed = toc(timer);
                [metrics,pointers] = feval(module,'metrics',plan);
                returned = feval(module,'pointer',destroyedUnique);
            case "inverse-destructive-aliased"
                destroyedAliased = complex(real(referenceSpectrum),imag(referenceSpectrum));
                aliasedInput = destroyedAliased;
                timer = tic;
                [destroyedAliased,inverseAliasedOutput] = feval(module,'inverseDestructive',plan,destroyedAliased,inverseAliasedOutput);
                elapsed = toc(timer);
                [metrics,pointers] = feval(module,'metrics',plan);
                returned = feval(module,'pointer',destroyedAliased);
            case "inverse-destructive-complex"
                destroyedComplex = complex(real(referenceSpectrum),imag(referenceSpectrum));
                wrapperOriginal = feval(module,'pointer',destroyedComplex);
                timer = tic;
                [destroyedComplex,inverseComplexOutput] = feval(module,'inverseDestructive',plan,complex(destroyedComplex),inverseComplexOutput);
                elapsed = toc(timer);
                [metrics,pointers] = feval(module,'metrics',plan);
                returned = feval(module,'pointer',destroyedComplex);
            case "chain-fftw-forward-destructive-inverse"
                timer = tic;
                chainSpectrum = feval(module,'forward',plan,x,'fftw-buffer');
                forwardElapsed = toc(timer);
                [forwardMetrics,forwardPointers] = feval(module,'metrics',plan);
                handoffToken = feval(module,'pointer',chainSpectrum);
                timer = tic;
                [chainSpectrum,chainOutput] = feval(module,'inverseDestructive',plan,chainSpectrum,chainOutput);
                inverseElapsed = toc(timer);
                [inverseMetrics,inversePointers] = feval(module,'metrics',plan);
                elapsed = forwardElapsed + inverseElapsed;
                metrics = combineMetrics(forwardMetrics,inverseMetrics);
                pointers = inversePointers;
                wrapperOriginal = handoffToken;
                returned = feval(module,'pointer',chainSpectrum);
        end
        if iRound > nWarmups
            sampleIndex = iRound-nWarmups;
            totalSamples(iOperation,sampleIndex) = elapsed;
            metricSamples(iOperation,:,sampleIndex) = metrics;
            pointerSamples(iOperation,:,sampleIndex) = pointers;
            returnedTokens(iOperation,sampleIndex) = returned;
            wrapperOriginalTokens(iOperation,sampleIndex) = wrapperOriginal;
        end
    end
end

noopSamples = nan(1,nSamples);
for iSample = 1:(nWarmups+nSamples)
    timer = tic;
    feval(module,'noop');
    elapsed = toc(timer);
    if iSample > nWarmups, noopSamples(iSample-nWarmups) = elapsed; end
end

actuals = {forwardFactory,forwardCaller,forwardComplex,forwardMatlabBuffer,forwardFFTWBuffer,scaleFactor*inverseAllocating,scaleFactor*inversePreallocated,scaleFactor*inverseUniqueOutput,scaleFactor*inverseAliasedOutput,scaleFactor*inverseComplexOutput,scaleFactor*chainOutput};
references = {referenceSpectrum,referenceSpectrum,referenceSpectrum,referenceSpectrum,referenceSpectrum,x,x,x,x,x,x};
operations = repmat(emptyOperation(),nOperations,1);
for iOperation = 1:nOperations
    [absoluteError,relativeError] = numericalError(actuals{iOperation},references{iOperation});
    operationMetricSamples = reshape(metricSamples(iOperation,:,:),17,nSamples).';
    operationPointerSamples = reshape(pointerSamples(iOperation,:,:),5,nSamples).';
    operations(iOperation) = operationRecord(operationIds(iOperation),totalSamples(iOperation,:),operationMetricSamples,operationPointerSamples,returnedTokens(iOperation,:),wrapperOriginalTokens(iOperation,:),absoluteError,relativeError,errorTolerance,prod(sz)*8,numel(referenceSpectrum)*16);
end

preservingUnchanged = isequaln(preservingInput,referenceSpectrum);
destructiveChanged = ~isequaln(destroyedUnique,referenceSpectrum);
aliasPreserved = isequaln(aliasedInput,referenceSpectrum);
clear actuals
clear forwardFFTWBuffer chainSpectrum
lifetime = feval(module,'lifetime');
lifetimeRecord.created = lifetime(1);
lifetimeRecord.freed = lifetime(2);
lifetimeRecord.outstanding = lifetime(3);
lifetimeRecord.balanced = lifetime(1) == lifetime(2) && lifetime(3) == 0;

workload = emptyWorkload();
workload.size = sz;
workload.seed = seed;
workload.configuration = configuration;
workload.complexSize = complexSize;
workload.scaleFactor = scaleFactor;
workload.plan.planningSeconds = planningSeconds;
workload.plan.inputAlignmentClass = inputAlignment;
workload.plan.outputAlignmentClass = outputAlignment;
workload.plan.memoryEstimate = planMemory;
workload.nWarmups = nWarmups;
workload.nSamples = nSamples;
workload.noopMexSamplesSeconds = noopSamples;
workload.noopMexMedianSeconds = median(noopSamples);
workload.operations = operations;
workload.correctness.maximumRelativeError = max([operations.maximumRelativeError]);
workload.correctness.passed = all([operations.correctnessPassed]);
workload.semantics.preservingInputUnchanged = preservingUnchanged;
workload.semantics.destructiveInputChanged = destructiveChanged;
workload.semantics.aliasPreserved = aliasPreserved;
workload.lifetime = lifetimeRecord;
workload.storage = storageRecord(sz,complexSize,planMemory);
workload.wrapperOverhead.forwardComplexMinusDirectSeconds = operations(3).medianSeconds-operations(2).medianSeconds;
workload.wrapperOverhead.inverseComplexMinusDirectSeconds = operations(10).medianSeconds-operations(8).medianSeconds;
if ~workload.correctness.passed || ~preservingUnchanged || ~destructiveChanged || ~aliasPreserved || ~lifetimeRecord.balanced
    error('FFTWMexOwnership:WorkloadValidationFailed','Validation failed for %s: correctness=%d, preserving=%d, destructive=%d, alias=%d, lifetime=%d (created=%g, freed=%g, outstanding=%g).',formatSize(sz),workload.correctness.passed,preservingUnchanged,destructiveChanged,aliasPreserved,lifetimeRecord.balanced,lifetimeRecord.created,lifetimeRecord.freed,lifetimeRecord.outstanding);
end
clear planCleanup
end

function memory = probePlanMemory(module,sz,transformOrder,threads,planner,alignmentMode,realTemplate,complexTemplate,nSamples,knownScratchBytes)
allocatorDeltas = nan(1,nSamples);
residentDeltas = nan(1,nSamples);
planningSamples = nan(1,nSamples);
for iSample = 1:nSamples
    before = feval(module,'memory');
    [temporaryPlan,~,~,planningSamples(iSample)] = feval(module,'create',sz,transformOrder,threads,planner,char(alignmentMode),realTemplate,complexTemplate);
    after = feval(module,'memory');
    allocatorDeltas(iSample) = max(0,after(1)-before(1)-knownScratchBytes);
    residentDeltas(iSample) = max(0,after(2)-before(2)-knownScratchBytes);
    feval(module,'free',temporaryPlan);
end
memory.method = "macOS default malloc-zone and resident-size deltas; known preserving scratch subtracted";
memory.knownPersistentScratchBytes = knownScratchBytes;
memory.allocatorDeltaSamplesBytes = allocatorDeltas;
memory.residentDeltaSamplesBytes = residentDeltas;
memory.planningSamplesSeconds = planningSamples;
memory.planOwnedAllocatorMedianBytes = median(allocatorDeltas,'omitnan');
memory.planOwnedAllocatorRangeBytes = [min(allocatorDeltas,[],'omitnan') max(allocatorDeltas,[],'omitnan')];
memory.planOwnedResidentMedianBytes = median(residentDeltas,'omitnan');
memory.planOwnedResidentRangeBytes = [min(residentDeltas,[],'omitnan') max(residentDeltas,[],'omitnan')];
memory.isEstimate = true;
end

function combined = combineMetrics(forward,inverse)
combined = forward + inverse;
combined(7:12) = forward(7:12) + inverse(7:12);
combined(13:15) = inverse(13:15);
combined(16) = max(forward(16),inverse(16));
combined(17) = inverse(17);
end

function operation = operationRecord(id,totalSamples,metricSamples,pointerSamples,returnedTokens,wrapperOriginalTokens,absoluteError,relativeError,errorTolerance,realBytes,halfSpectrumBytes)
operation = emptyOperation();
operation.id = id;
operation.direction = operationDirection(id);
operation.totalMexSamplesSeconds = totalSamples;
operation.medianSeconds = median(totalSamples);
metricNames = ["allocation","bufferWrap","explicitMemcpy","mutableDetach","kernel","internalPipeline"];
for iMetric = 1:numel(metricNames)
    field = metricNames(iMetric) + "SamplesSeconds";
    operation.(field) = metricSamples(:,iMetric).';
    medianField = metricNames(iMetric) + "MedianSeconds";
    operation.(medianField) = median(metricSamples(:,iMetric));
end
operation.boundaryResidualSamplesSeconds = max(0,totalSamples-metricSamples(:,6).');
operation.boundaryResidualMedianSeconds = median(operation.boundaryResidualSamplesSeconds);
operation.allocationCountSamples = metricSamples(:,7).';
operation.allocatedBytesSamples = metricSamples(:,8).';
operation.explicitCopyCountSamples = metricSamples(:,9).';
operation.explicitCopiedBytesSamples = metricSamples(:,10).';
operation.detectedCopyCountSamples = metricSamples(:,11).';
operation.detectedCopiedBytesSamples = metricSamples(:,12).';
operation.inputAlignmentClasses = metricSamples(:,13).';
operation.outputAlignmentClasses = metricSamples(:,14).';
operation.destroyedInputSamples = logical(metricSamples(:,15).');
operation.persistentScratchBytes = max(metricSamples(:,16));
operation.outstandingFftwBuffersAtReturn = metricSamples(:,17).';
operation.inputBeforeTokens = pointerSamples(:,1).';
operation.inputMutableTokens = pointerSamples(:,2).';
operation.outputBeforeTokens = pointerSamples(:,3).';
operation.outputMutableTokens = pointerSamples(:,4).';
operation.wrappedTokens = pointerSamples(:,5).';
operation.returnedTokens = returnedTokens;
operation.wrapperOriginalTokens = wrapperOriginalTokens;
if startsWith(id,"inverse-destructive-") || startsWith(id,"chain-")
    operation.returnedPointerPreserved = all(returnedTokens == operation.inputMutableTokens);
else
    operation.returnedPointerPreserved = all(returnedTokens == operation.wrappedTokens);
end
operation.inputDetached = operation.inputBeforeTokens ~= operation.inputMutableTokens;
operation.outputDetached = operation.outputBeforeTokens ~= operation.outputMutableTokens & operation.outputBeforeTokens ~= 0;
operation.wrapperAllocated = wrapperOriginalTokens ~= 0 & wrapperOriginalTokens ~= operation.outputBeforeTokens & wrapperOriginalTokens ~= operation.inputBeforeTokens;
if startsWith(id,"chain-")
    operation.handoffPointerPreserved = all(wrapperOriginalTokens == operation.inputBeforeTokens);
else
    operation.handoffPointerPreserved = [];
end
operation.maximumAbsoluteError = absoluteError;
operation.maximumRelativeError = relativeError;
operation.correctnessPassed = relativeError <= errorTolerance;
[operation.allocationOwner,operation.callerPreallocatedBytes] = allocationLedger(id,realBytes,halfSpectrumBytes);
end

function [owner,callerBytes] = allocationLedger(id,realBytes,halfSpectrumBytes)
callerBytes = 0;
switch id
    case "forward-factory-array"
        owner = "MATLAB Data API TypedArray";
    case {"forward-caller-direct","forward-caller-complex"}
        owner = "caller; MATLAB copy-on-write only if detached";
        callerBytes = halfSpectrumBytes;
    case "forward-matlab-buffer"
        owner = "MATLAB Data API buffer";
    case "forward-fftw-buffer"
        owner = "FFTW buffer with custom deleter";
    case "inverse-preserving-allocating"
        owner = "MATLAB Data API real output plus persistent MEX scratch";
    case "inverse-preserving-preallocated"
        owner = "caller real output plus persistent MEX scratch";
        callerBytes = realBytes;
    case {"inverse-destructive-unique","inverse-destructive-aliased","inverse-destructive-complex"}
        owner = "caller spectrum and real output; MATLAB copy-on-write if detached";
        callerBytes = halfSpectrumBytes + realBytes;
    case "chain-fftw-forward-destructive-inverse"
        owner = "FFTW spectrum with custom deleter plus caller real output";
        callerBytes = realBytes;
end
end

function direction = operationDirection(id)
if startsWith(id,"forward-")
    direction = "forward";
elseif startsWith(id,"inverse-")
    direction = "inverse";
else
    direction = "roundtrip";
end
end

function storage = storageRecord(sz,complexSize,planMemory)
realElements = prod(sz);
halfElements = prod(complexSize);
halfXElements = (floor(sz(1)/2)+1)*sz(2)*sz(3);
halfYElements = sz(1)*(floor(sz(2)/2)+1)*sz(3);
storage.realArrayBytes = realElements*8;
storage.matlabFullSpectrumBytes = realElements*16;
storage.halfSpectrumBytes = halfElements*16;
storage.halfToFullSpectrumRatio = storage.halfSpectrumBytes/storage.matlabFullSpectrumBytes;
storage.standaloneLegacy.wrapperScratchBytes = storage.halfSpectrumBytes;
storage.standaloneLegacy.mexPreservingScratchBytes = storage.halfSpectrumBytes;
storage.standaloneLegacy.totalKnownBytes = 2*storage.halfSpectrumBytes;
storage.waveVortexMatlab.fullComplexBufferBytes = realElements*16;
storage.waveVortexMatlab.totalKnownBytes = storage.waveVortexMatlab.fullComplexBufferBytes;
storage.waveVortexFFTW.fullComplexBufferBytes = realElements*16;
storage.waveVortexFFTW.realBufferBytes = realElements*8;
storage.waveVortexFFTW.explicitHalfSpectrumBuffersBytes = 16*(halfXElements+2*halfYElements);
storage.waveVortexFFTW.wrapperScratchBuffersBytes = storage.waveVortexFFTW.explicitHalfSpectrumBuffersBytes;
storage.waveVortexFFTW.mexPreservingBuffersBytes = storage.waveVortexFFTW.explicitHalfSpectrumBuffersBytes;
storage.waveVortexFFTW.totalKnownBytes = storage.waveVortexFFTW.fullComplexBufferBytes + storage.waveVortexFFTW.realBufferBytes + 3*storage.waveVortexFFTW.explicitHalfSpectrumBuffersBytes;
storage.leanBackend.destructiveOnlyPersistentSpectrumBytes = 0;
storage.leanBackend.withPreservingInverseScratchBytes = storage.halfSpectrumBytes;
storage.estimatedPlanOwnedAllocatorBytes = planMemory.planOwnedAllocatorMedianBytes;
end

function recommendation = chooseRecommendation(workloads,gateSizes,errorTolerance)
forwardIds = ["forward-factory-array","forward-caller-direct","forward-caller-complex","forward-matlab-buffer","forward-fftw-buffer"];
nModes = numel(forwardIds);
nGates = size(gateSizes,1);
times = inf(nModes,nGates);
normalizedTimes = inf(nModes,nGates);
valid = true(nModes,nGates);
for iGate = 1:nGates
    workloadIndex = find(arrayfun(@(workload) isequal(double(workload.size),gateSizes(iGate,:)),workloads),1);
    if isempty(workloadIndex), error('FFTWMexOwnership:GateWorkloadMissing','A required gate workload was not measured.'); end
    gateOperations = repmat(emptyOperation(),nModes,1);
    for iMode = 1:nModes
        operation = workloads(workloadIndex).operations([workloads(workloadIndex).operations.id] == forwardIds(iMode));
        gateOperations(iMode) = operation;
        times(iMode,iGate) = operation.medianSeconds;
        valid(iMode,iGate) = operation.maximumRelativeError <= errorTolerance && operation.returnedPointerPreserved && all(operation.detectedCopiedBytesSamples == 0) && workloads(workloadIndex).lifetime.balanced;
    end
    commonKernelMedian = median([gateOperations.kernelSamplesSeconds]);
    for iMode = 1:nModes
        normalizedTimes(iMode,iGate) = commonKernelMedian + median(gateOperations(iMode).totalMexSamplesSeconds-gateOperations(iMode).kernelSamplesSeconds);
    end
end
gateMinimum = min(times,[],1);
eligible = all(valid,2) & all(times <= 1.05*gateMinimum,2);
selectionTimes = times;
selectionBasis = "raw complete MEX medians";
strictRawFivePercentCandidateExists = any(eligible);
if ~strictRawFivePercentCandidateExists
    normalizedMinimum = min(normalizedTimes,[],1);
    eligible = all(valid,2) & all(normalizedTimes <= 1.05*normalizedMinimum,2);
    selectionTimes = normalizedTimes;
    selectionBasis = "common FFT-kernel median plus measured ownership overhead because raw kernel variance left no common 5 percent candidate";
end
geometricMeans = exp(mean(log(selectionTimes),2));
geometricMeans(~eligible) = Inf;
best = min(geometricMeans);
nearBest = find(eligible & geometricMeans <= 1.03*best);
preference = [2 1 4 3 5];
selected = [];
for candidate = preference
    if ismember(candidate,nearBest)
        selected = candidate;
        break
    end
end
if isempty(selected), error('FFTWMexOwnership:NoValidForwardOwnership','No forward ownership model satisfied the gate selection rule.'); end

inverseValid = true(1,nGates);
for iGate = 1:nGates
    workloadIndex = find(arrayfun(@(workload) isequal(double(workload.size),gateSizes(iGate,:)),workloads),1);
    operation = workloads(workloadIndex).operations([workloads(workloadIndex).operations.id] == "inverse-destructive-unique");
    inverseValid(iGate) = operation.maximumRelativeError <= errorTolerance && ~any(operation.inputDetached) && all(operation.explicitCopiedBytesSamples == 0);
end
recommendation.forwardOwnership = forwardIds(selected);
recommendation.forwardGateTimesSeconds = times(selected,:);
recommendation.forwardSelectionGateTimesSeconds = selectionTimes(selected,:);
recommendation.forwardGeometricMeanSeconds = geometricMeans(selected);
recommendation.forwardSelectionBasis = selectionBasis;
recommendation.strictRawFivePercentCandidateExists = strictRawFivePercentCandidateExists;
recommendation.forwardRationale = "Valid zero-copy model selected under the 5 percent gate and 3 percent ownership-simplicity preference rules.";
switch recommendation.forwardOwnership
    case "forward-caller-direct"
        recommendation.forwardContract = "Reuse one caller-owned half-spectrum buffer, move it directly into MEX, reassign the returned array, and omit complex(...).";
    case "forward-factory-array"
        recommendation.forwardContract = "Allocate the half-spectrum as an ordinary MATLAB Data API TypedArray inside MEX.";
    case "forward-matlab-buffer"
        recommendation.forwardContract = "Allocate a MATLAB-managed buffer inside MEX and return it through createArrayFromBuffer.";
    case "forward-caller-complex"
        recommendation.forwardContract = "Reuse a caller-owned half-spectrum buffer through complex(...); the measured wrapper must remain allocation-free.";
    case "forward-fftw-buffer"
        recommendation.forwardContract = "Allocate with FFTW and transfer ownership through createArrayFromBuffer with the fftw_free custom deleter.";
end
recommendation.destructiveInverse = all(inverseValid);
if recommendation.destructiveInverse
    recommendation.inverseContract = "Move a uniquely owned half-spectrum directly, return and reassign the destroyed spectrum, and omit complex(...).";
else
    recommendation.inverseContract = "No copy-free destructive inverse contract was demonstrated on both gate workloads.";
end
recommendation.aliasContract = "Aliased spectra may detach through MATLAB copy-on-write; aliasing is outside the zero-copy contract.";
recommendation.preservingInverseContract = "Input-preserving c2r retains exactly one explicit half-spectrum memcpy into persistent scratch.";
recommendation.doesNotDeclareMilestoneDecision = true;
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
end
end

function output = fftAlong(input,order)
output = input;
for dimension = order, output = fft(output,[],dimension); end
end

function [maximumAbsoluteError,maximumRelativeError] = numericalError(actual,reference)
maximumAbsoluteError = max(abs(actual-reference),[],'all');
referenceMagnitude = max(abs(reference),[],'all');
if referenceMagnitude == 0
    maximumRelativeError = maximumAbsoluteError;
else
    maximumRelativeError = maximumAbsoluteError/referenceMagnitude;
end
end

function writeResultArtifacts(result,runDirectory)
writeText(fullfile(runDirectory,"ownership-benchmark.json"),jsonencode(result));
writeText(fullfile(runDirectory,"summary.md"),markdownSummary(result));
end

function summary = markdownSummary(result)
lines = strings(0,1);
lines(end+1) = "# FFTW MEX ownership, alignment, and copy benchmark";
lines(end+1) = "";
lines(end+1) = "Status: **" + upper(result.status) + "**";
lines(end+1) = "";
lines(end+1) = "This issue #39 report characterizes ownership and copies without making the milestone's GO or NO-GO decision.";
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
lines(end+1) = markdownRow("Physical memory",sprintf('%.3f GB',result.environment.physicalMemoryBytes/1024^3));
if result.status == "passed"
    lines(end+1) = "";
    lines(end+1) = "## Ownership timing";
    lines(end+1) = "";
    lines(end+1) = "Times are medians in milliseconds. Boundary residual is complete MEX time minus the instrumented internal pipeline.";
    lines(end+1) = "";
    lines(end+1) = "| Size | Operation | Total | Kernel | Allocate | Wrap | memcpy | Detach | Boundary | Detected copied MiB |";
    lines(end+1) = "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|";
    for iWorkload = 1:numel(result.workloads)
        for iOperation = 1:numel(result.workloads(iWorkload).operations)
            operation = result.workloads(iWorkload).operations(iOperation);
            copiedMiB = median(operation.detectedCopiedBytesSamples)/1024^2;
            lines(end+1) = sprintf('| %s | %s | %.4f | %.4f | %.4f | %.4f | %.4f | %.4f | %.4f | %.3f |',formatSize(result.workloads(iWorkload).size),operation.id,1e3*operation.medianSeconds,1e3*operation.kernelMedianSeconds,1e3*operation.allocationMedianSeconds,1e3*operation.bufferWrapMedianSeconds,1e3*operation.explicitMemcpyMedianSeconds,1e3*operation.mutableDetachMedianSeconds,1e3*operation.boundaryResidualMedianSeconds,copiedMiB); %#ok<AGROW>
        end
    end
    lines(end+1) = "";
    lines(end+1) = "### Allocation and copy ledger";
    lines(end+1) = "";
    lines(end+1) = "| Size | Operation | Allocation owner | Timed allocated MiB | Explicit copied MiB | Detected COW MiB | Caller preallocated MiB | Persistent scratch MiB |";
    lines(end+1) = "|---|---|---|---:|---:|---:|---:|---:|";
    for iWorkload = 1:numel(result.workloads)
        for iOperation = 1:numel(result.workloads(iWorkload).operations)
            operation = result.workloads(iWorkload).operations(iOperation);
            lines(end+1) = sprintf('| %s | %s | %s | %.3f | %.3f | %.3f | %.3f | %.3f |',formatSize(result.workloads(iWorkload).size),operation.id,operation.allocationOwner,toMiB(median(operation.allocatedBytesSamples)),toMiB(median(operation.explicitCopiedBytesSamples)),toMiB(median(operation.detectedCopiedBytesSamples)),toMiB(operation.callerPreallocatedBytes),toMiB(operation.persistentScratchBytes)); %#ok<AGROW>
        end
    end
    lines(end+1) = "";
    lines(end+1) = "## Zero-copy and lifetime evidence";
    lines(end+1) = "";
    lines(end+1) = "| Size | Operation | Input detached | Output detached | Wrapper allocated | Return pointer preserved | Explicit copied MiB | Error |";
    lines(end+1) = "|---|---|:---:|:---:|:---:|:---:|---:|---:|";
    for iWorkload = 1:numel(result.workloads)
        for iOperation = 1:numel(result.workloads(iWorkload).operations)
            operation = result.workloads(iWorkload).operations(iOperation);
            lines(end+1) = sprintf('| %s | %s | %s | %s | %s | %s | %.3f | %.3g |',formatSize(result.workloads(iWorkload).size),operation.id,yesNo(any(operation.inputDetached)),yesNo(any(operation.outputDetached)),yesNo(any(operation.wrapperAllocated)),yesNo(operation.returnedPointerPreserved),median(operation.explicitCopiedBytesSamples)/1024^2,operation.maximumRelativeError); %#ok<AGROW>
        end
    end
    lines(end+1) = "";
    lines(end+1) = "| Size | FFTW buffers created | Freed | Outstanding | Balanced |";
    lines(end+1) = "|---|---:|---:|---:|:---:|";
    for iWorkload = 1:numel(result.workloads)
        lifetime = result.workloads(iWorkload).lifetime;
        lines(end+1) = sprintf('| %s | %g | %g | %g | %s |',formatSize(result.workloads(iWorkload).size),lifetime.created,lifetime.freed,lifetime.outstanding,yesNo(lifetime.balanced)); %#ok<AGROW>
    end
    lines(end+1) = "";
    lines(end+1) = "## Persistent storage";
    lines(end+1) = "";
    lines(end+1) = "| Size | MATLAB full spectrum MiB | Half spectrum MiB | Legacy standalone MiB | WaveVortex MATLAB MiB | Current FFTW backend MiB | Lean preserving scratch MiB | Estimated plan MiB |";
    lines(end+1) = "|---|---:|---:|---:|---:|---:|---:|---:|";
    for iWorkload = 1:numel(result.workloads)
        storage = result.workloads(iWorkload).storage;
        lines(end+1) = sprintf('| %s | %.3f | %.3f | %.3f | %.3f | %.3f | %.3f | %.3f |',formatSize(result.workloads(iWorkload).size),toMiB(storage.matlabFullSpectrumBytes),toMiB(storage.halfSpectrumBytes),toMiB(storage.standaloneLegacy.totalKnownBytes),toMiB(storage.waveVortexMatlab.totalKnownBytes),toMiB(storage.waveVortexFFTW.totalKnownBytes),toMiB(storage.leanBackend.withPreservingInverseScratchBytes),toMiB(storage.estimatedPlanOwnedAllocatorBytes)); %#ok<AGROW>
    end
    lines(end+1) = "";
    lines(end+1) = "Plan-owned values are allocator-delta estimates; application-owned buffers and copies are exact.";
    lines(end+1) = "";
    lines(end+1) = "## Recommendation";
    lines(end+1) = "";
    lines(end+1) = "- Forward ownership: `" + result.recommendation.forwardOwnership + "`.";
    lines(end+1) = "- Forward contract: " + result.recommendation.forwardContract;
    lines(end+1) = "- Selection basis: " + result.recommendation.forwardSelectionBasis + ".";
    lines(end+1) = "- A common raw 5% candidate existed: **" + yesNo(result.recommendation.strictRawFivePercentCandidateExists) + "**.";
    lines(end+1) = "- Destructive inverse demonstrated without a spectrum copy: **" + yesNo(result.recommendation.destructiveInverse) + "**.";
    lines(end+1) = "- Contract: " + result.recommendation.inverseContract;
    lines(end+1) = "- Preserving inverse: " + result.recommendation.preservingInverseContract;
end
if result.status == "failed" && ~isempty(result.failure)
    lines(end+1) = "";
    lines(end+1) = "## Failure";
    lines(end+1) = "";
    lines(end+1) = "- Stage: `" + result.failure.stage + "`";
    lines(end+1) = "- Identifier: `" + result.failure.identifier + "`";
    lines(end+1) = "- Message: " + result.failure.message;
end
summary = strjoin(lines,newline);
end

function writeText(path,text)
[fileId,message] = fopen(path,'w');
if fileId < 0, error('FFTWMexOwnership:ArtifactOpenFailed','Unable to open %s: %s',path,message); end
cleanup = onCleanup(@() fclose(fileId));
fwrite(fileId,char(text),'char');
clear cleanup
end

function environment = collectEnvironment()
environment.matlabVersion = string(version);
environment.matlabRelease = "R" + string(version('-release'));
environment.architecture = string(computer('arch'));
[~,operatingSystem] = system('uname -srvmp');
environment.operatingSystem = strtrim(string(operatingSystem));
[~,processor] = system('sysctl -n machdep.cpu.brand_string');
environment.processor = strtrim(string(processor));
[~,memory] = system('sysctl -n hw.memsize');
environment.physicalMemoryBytes = str2double(strtrim(memory));
compiler = mex.getCompilerConfigurations('C++','Selected');
if isempty(compiler)
    environment.compiler = "unknown";
else
    environment.compiler = string(compiler.Name) + " " + string(compiler.Version);
end
environment.hardwareThreads = maxNumCompThreads;
end

function sources = collectSourceProvenance(sourceDirectory)
names = ["runFFTWMexOwnershipBenchmark.m","buildFFTWMexOwnershipBenchmark.m","fftw_ownership_benchmark.cpp","runFFTWEngineLayoutBenchmark.m","fftw_engine_benchmark.cpp"];
sources = repmat(struct('path',"",'sha256',""),numel(names),1);
for iSource = 1:numel(names)
    path = fullfile(sourceDirectory,names(iSource));
    sources(iSource).path = names(iSource);
    sources(iSource).sha256 = fileSHA256(path);
end
end

function hash = fileSHA256(path)
bytes = fileread(path);
digest = java.security.MessageDigest.getInstance('SHA-256');
digest.update(uint8(bytes));
hashBytes = typecast(digest.digest(),'uint8');
hash = string(lower(reshape(dec2hex(hashBytes,2).',1,[])));
end

function restoreMatlabState(planner,wisdom,threads,randomState)
try
    fftw('dwisdom',[]);
    fftw('dwisdom',wisdom);
catch
end
try
    fftw('planner',char(planner));
catch
end
try
    maxNumCompThreads(threads);
catch
end
try
    rng(randomState);
catch
end
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

function row = markdownRow(field,value)
row = "| " + string(field) + " | " + replace(replace(string(value),"|","\|"),newline,"<br>") + " |";
end

function value = formatSize(sz)
value = strjoin(string(sz)," x ");
end

function value = yesNo(condition)
if condition, value = "yes"; else, value = "no"; end
end

function value = toMiB(bytes)
value = bytes/1024^2;
end

function engine = emptyEngine()
engine = struct('id',"",'module',"",'status',"",'version',"",'library',"",'alignmentMatchedAccepted',false,'alignmentMismatchRejected',false,'unalignedAccepted',false);
end

function configuration = emptyConfiguration()
configuration = struct('size',[],'engine',"",'module',"",'transformOrder',[],'layout',"",'strategy',"",'planner',"",'threads',[],'alignmentMode',"");
end

function operation = emptyOperation()
operation = struct('id',"",'direction',"",'totalMexSamplesSeconds',[],'medianSeconds',NaN,'allocationSamplesSeconds',[],'allocationMedianSeconds',NaN,'bufferWrapSamplesSeconds',[],'bufferWrapMedianSeconds',NaN,'explicitMemcpySamplesSeconds',[],'explicitMemcpyMedianSeconds',NaN,'mutableDetachSamplesSeconds',[],'mutableDetachMedianSeconds',NaN,'kernelSamplesSeconds',[],'kernelMedianSeconds',NaN,'internalPipelineSamplesSeconds',[],'internalPipelineMedianSeconds',NaN,'boundaryResidualSamplesSeconds',[],'boundaryResidualMedianSeconds',NaN,'allocationCountSamples',[],'allocatedBytesSamples',[],'explicitCopyCountSamples',[],'explicitCopiedBytesSamples',[],'detectedCopyCountSamples',[],'detectedCopiedBytesSamples',[],'inputAlignmentClasses',[],'outputAlignmentClasses',[],'destroyedInputSamples',[],'persistentScratchBytes',0,'outstandingFftwBuffersAtReturn',[],'inputBeforeTokens',uint64([]),'inputMutableTokens',uint64([]),'outputBeforeTokens',uint64([]),'outputMutableTokens',uint64([]),'wrappedTokens',uint64([]),'returnedTokens',uint64([]),'wrapperOriginalTokens',uint64([]),'returnedPointerPreserved',false,'inputDetached',[],'outputDetached',[],'wrapperAllocated',[],'handoffPointerPreserved',[],'allocationOwner',"",'callerPreallocatedBytes',0,'maximumAbsoluteError',NaN,'maximumRelativeError',NaN,'correctnessPassed',false);
end

function workload = emptyWorkload()
workload = struct('size',[],'seed',[],'configuration',emptyConfiguration(),'complexSize',[],'scaleFactor',NaN,'plan',struct,'nWarmups',[],'nSamples',[],'noopMexSamplesSeconds',[],'noopMexMedianSeconds',NaN,'operations',repmat(emptyOperation(),0,1),'correctness',struct,'semantics',struct,'lifetime',struct,'storage',struct,'wrapperOverhead',struct);
end

function recommendation = emptyRecommendation()
recommendation = struct('forwardOwnership',"",'forwardGateTimesSeconds',[],'forwardSelectionGateTimesSeconds',[],'forwardGeometricMeanSeconds',NaN,'forwardSelectionBasis',"",'strictRawFivePercentCandidateExists',false,'forwardRationale',"",'forwardContract',"",'destructiveInverse',false,'inverseContract',"",'aliasContract',"",'preservingInverseContract',"",'doesNotDeclareMilestoneDecision',true);
end
