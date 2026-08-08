function result = runFFTWR2RBenchmark(options)
% Benchmark normalized FFTW DCT-I/DST-I transforms and derive eligibility.
arguments
    options.outputDirectory (1,1) string = fullfile(fileparts(mfilename('fullpath')),"results","issue43")
    options.Nz (1,:) double = [33 65 129 257 513]
    options.batchCounts (1,:) double = [1 8320 33024 131584]
    options.dataTypes (1,:) string = ["real","complex"]
    options.transformTypes (1,:) string = ["cosine","sine"]
    options.planner (1,1) string = "measure"
    options.threadCount (1,1) double = maxNumCompThreads
    options.alignmentMode (1,1) string = "unaligned"
    options.plannerTimeLimitSeconds (1,1) double = 10
    options.nWarmups (1,1) double = 2
    options.nSamples (1,1) double = 7
    options.nSamplesLargest (1,1) double = 3
    options.seed (1,1) double = 43
    options.errorTolerance (1,1) double = 1e-12
    options.speedThreshold (1,1) double = 1.10
    options.shouldBuild (1,1) logical = true
    options.shouldWriteArtifacts (1,1) logical = true
    options.requireCanonicalPlatform (1,1) logical = true
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
fftw('planner',char(options.planner));
fftw('dwisdom',[]);

result = initializeResult(runId,runDirectory,sourceDirectory,options);
activeStage = "build";
try
    if options.shouldBuild
        RealToRealTransform.makeMexFiles;
        result.build.status = "built";
    elseif exist('fftw_r2r','file') ~= 3
        error('FFTWR2RBenchmark:MexMissing','fftw_r2r is unavailable and shouldBuild is false.');
    else
        result.build.status = "existing";
    end
    [versionText,libraryPath] = fftw_r2r('info');
    result.engine.version = string(versionText);
    result.engine.library = string(libraryPath);
    result.engine.module = "fftw_r2r";

    activeStage = "measurements";
    iConfiguration = 0;
    for n = options.Nz
        for batchCount = options.batchCounts
            for dataType = options.dataTypes
                for transformType = options.transformTypes
                    iConfiguration = iConfiguration+1;
                    nSamples = options.nSamples;
                    if n == options.Nz(end) && batchCount == options.batchCounts(end)
                        nSamples = options.nSamplesLargest;
                    end
                    rng(options.seed+iConfiguration-1,'twister');
                    records = benchmarkConfiguration(n,batchCount,dataType,transformType,options,nSamples);
                    result.workloads = [result.workloads; records];
                end
            end
        end
    end

    activeStage = "eligibility";
    result.eligibility = deriveEligibility(result.workloads,options);
    result.status = "passed";
    result.completedAtUTC = utcTimestamp;
    if options.shouldWriteArtifacts
        writeResultArtifacts(result,runDirectory);
    end
catch exception
    result.status = "failed";
    result.completedAtUTC = utcTimestamp;
    result.failure = failureRecord(exception,activeStage);
    if options.shouldWriteArtifacts
        try
            writeResultArtifacts(result,runDirectory);
        catch artifactException
            warning('FFTWR2RBenchmark:FailureArtifactWriteFailed','Unable to write issue #43 failure artifacts: %s',artifactException.message);
        end
    end
    rethrow(exception);
end
clear stateCleanup
end

function validateOptions(options)
validateattributes(options.Nz,{'double'},{'row','integer','>=',3},mfilename,'Nz');
validateattributes(options.batchCounts,{'double'},{'row','integer','positive'},mfilename,'batchCounts');
validateattributes(options.threadCount,{'double'},{'scalar','integer','positive'},mfilename,'threadCount');
validateattributes(options.plannerTimeLimitSeconds,{'double'},{'scalar','positive','finite'},mfilename,'plannerTimeLimitSeconds');
validateattributes(options.nWarmups,{'double'},{'scalar','integer','nonnegative'},mfilename,'nWarmups');
validateattributes(options.nSamples,{'double'},{'scalar','integer','positive'},mfilename,'nSamples');
validateattributes(options.nSamplesLargest,{'double'},{'scalar','integer','positive'},mfilename,'nSamplesLargest');
validateattributes(options.seed,{'double'},{'scalar','integer','nonnegative'},mfilename,'seed');
validateattributes(options.errorTolerance,{'double'},{'scalar','nonnegative','finite'},mfilename,'errorTolerance');
validateattributes(options.speedThreshold,{'double'},{'scalar','>',1,'finite'},mfilename,'speedThreshold');
if any(~ismember(options.dataTypes,["real","complex"]))
    error('FFTWR2RBenchmark:InvalidDataType','dataTypes must contain real or complex.');
end
if any(~ismember(options.transformTypes,["cosine","sine"]))
    error('FFTWR2RBenchmark:InvalidTransformType','transformTypes must contain cosine or sine.');
end
if ~ismember(options.planner,["estimate","measure","patient","exhaustive"])
    error('FFTWR2RBenchmark:InvalidPlanner','Unknown FFTW planner.');
end
if ~ismember(options.alignmentMode,["matched","unaligned"])
    error('FFTWR2RBenchmark:InvalidAlignmentMode','Unknown alignment mode.');
end
if options.requireCanonicalPlatform && (~startsWith(string(version('-release')),"2026a",IgnoreCase=true) || string(computer('arch')) ~= "maca64")
    error('FFTWR2RBenchmark:UnsupportedPlatform','The canonical issue #43 benchmark targets MATLAB R2026a on maca64.');
end
end

function result = initializeResult(runId,runDirectory,sourceDirectory,options)
result.schemaVersion = "1.0.0";
result.status = "running";
result.runId = runId;
result.generatedAtUTC = utcTimestamp;
result.completedAtUTC = "";
result.environment = collectEnvironment;
result.configuration.Nz = options.Nz;
result.configuration.batchCounts = options.batchCounts;
result.configuration.batchLabels = batchLabels(options.batchCounts);
result.configuration.dataTypes = options.dataTypes;
result.configuration.transformTypes = options.transformTypes;
result.configuration.directions = ["forward","inverse"];
result.configuration.planner = options.planner;
result.configuration.requestedThreads = options.threadCount;
result.configuration.observedThreads = min(options.threadCount,maxNumCompThreads);
result.configuration.alignmentMode = options.alignmentMode;
result.configuration.plannerTimeLimitSeconds = options.plannerTimeLimitSeconds;
result.configuration.nWarmups = options.nWarmups;
result.configuration.nSamples = options.nSamples;
result.configuration.nSamplesLargest = options.nSamplesLargest;
result.configuration.seed = options.seed;
result.configuration.errorTolerance = options.errorTolerance;
result.configuration.speedThreshold = options.speedThreshold;
result.configuration.operationOrder = "four-method round-robin with rotating first operation";
result.configuration.eligibilityPolicy = "bounded intervals joining adjacent passing batch anchors; no extrapolation";
result.sources = collectSourceProvenance(sourceDirectory);
result.build = struct;
result.engine = struct('module',"",'version',"",'library',"");
result.workloads = repmat(emptyWorkload,0,1);
result.eligibility = repmat(emptyEligibility,0,1);
result.failure = [];
result.artifacts.directory = runDirectory;
result.artifacts.json = "real-to-real-benchmark.json";
result.artifacts.markdown = "summary.md";
end

function records = benchmarkConfiguration(n,batchCount,dataType,transformType,options,nSamples)
sz = [n batchCount];
if dataType == "real"
    physical = randn(sz);
else
    physical = complex(randn(sz),randn(sz));
end
if transformType == "sine"
    physical([1 end],:) = 0;
end
[forwardMatrix,backMatrix] = referenceMatrices(n,transformType);
coefficients = forwardMatrix*physical;
transform = RealToRealTransform(sz,dims=1,transform=transformType,dataType=dataType,planner=options.planner,nCores=min(options.threadCount,maxNumCompThreads),alignmentMode=options.alignmentMode,plannerTimeLimitSeconds=options.plannerTimeLimitSeconds);
transformCleanup = onCleanup(@() delete(transform));

forwardSetup.input = physical;
forwardSetup.reference = coefficients;
forwardSetup.preallocated = zerosLike(transform.spectralSize,dataType);
forwardSetup.coordinate = (0:n-1)'/(n-1);
forwardSetup.frequency = [];
inverseSetup.input = coefficients;
inverseSetup.reference = backMatrix*coefficients;
inverseSetup.preallocated = zerosLike(transform.realSize,dataType);
if transformType == "cosine"
    inverseSetup.frequency = (0:n-1)'/(2*(n-1));
else
    inverseSetup.frequency = (1:n-2)'/(2*(n-1));
end
inverseSetup.coordinate = [];

records = repmat(emptyWorkload,2,1);
records(1) = benchmarkDirection(transform,forwardMatrix,backMatrix,forwardSetup,"forward",n,batchCount,dataType,transformType,options,nSamples);
records(2) = benchmarkDirection(transform,forwardMatrix,backMatrix,inverseSetup,"inverse",n,batchCount,dataType,transformType,options,nSamples);
clear transformCleanup
end

function workload = benchmarkDirection(transform,forwardMatrix,backMatrix,setup,direction,n,batchCount,dataType,transformType,options,nSamples)
methodIds = ["dense-matrix","fft-extension","fftw-allocating","fftw-preallocated"];
nMethods = numel(methodIds);
nRounds = options.nWarmups+nSamples;
totalSamples = nan(nMethods,nSamples);
metricSamples = nan(nMethods,13,nSamples);
outputs = cell(nMethods,1);
preallocated = setup.preallocated;
for iRound = 1:nRounds
    firstMethod = mod(iRound-1,nMethods)+1;
    methodOrder = [firstMethod:nMethods 1:firstMethod-1];
    for iMethod = methodOrder
        switch methodIds(iMethod)
            case "dense-matrix"
                timer = tic;
                if direction == "forward"
                    outputs{iMethod} = forwardMatrix*setup.input;
                else
                    outputs{iMethod} = backMatrix*setup.input;
                end
                elapsed = toc(timer);
                metrics = nan(1,13);
            case "fft-extension"
                timer = tic;
                outputs{iMethod} = extensionTransform(setup.input,setup.coordinate,setup.frequency,transformType,direction);
                elapsed = toc(timer);
                metrics = nan(1,13);
            case "fftw-allocating"
                timer = tic;
                if direction == "forward"
                    outputs{iMethod} = transform.transformForward(setup.input);
                else
                    outputs{iMethod} = transform.transformBack(setup.input);
                end
                elapsed = toc(timer);
                metrics = [];
            case "fftw-preallocated"
                timer = tic;
                if direction == "forward"
                    preallocated = transform.transformForwardIntoArray(setup.input,preallocated);
                else
                    preallocated = transform.transformBackIntoArray(setup.input,preallocated);
                end
                elapsed = toc(timer);
                outputs{iMethod} = preallocated;
                metrics = [];
        end
        if startsWith(methodIds(iMethod),"fftw")
            metrics = transformMetrics(transform);
        end
        if iRound > options.nWarmups
            sampleIndex = iRound-options.nWarmups;
            totalSamples(iMethod,sampleIndex) = elapsed;
            metricSamples(iMethod,:,sampleIndex) = metrics;
        end
    end
end

operations = repmat(emptyOperation,nMethods,1);
for iMethod = 1:nMethods
    operations(iMethod).id = methodIds(iMethod);
    operations(iMethod).totalSamplesSeconds = totalSamples(iMethod,:);
    operations(iMethod).medianSeconds = median(totalSamples(iMethod,:));
    if startsWith(methodIds(iMethod),"fftw")
        values = reshape(metricSamples(iMethod,:,:),13,nSamples).';
        operations(iMethod).allocationSamplesSeconds = values(:,1).';
        operations(iMethod).wrapSamplesSeconds = values(:,2).';
        operations(iMethod).kernelSamplesSeconds = values(:,3).';
        operations(iMethod).normalizationSamplesSeconds = values(:,4).';
        operations(iMethod).detachSamplesSeconds = values(:,5).';
        operations(iMethod).internalSamplesSeconds = values(:,6).';
        operations(iMethod).pipelineSamplesSeconds = (values(:,3)+values(:,4)).';
        operations(iMethod).boundaryResidualSamplesSeconds = max(0,totalSamples(iMethod,:)-values(:,6).');
        operations(iMethod).allocationMedianSeconds = median(values(:,1));
        operations(iMethod).wrapMedianSeconds = median(values(:,2));
        operations(iMethod).kernelMedianSeconds = median(values(:,3));
        operations(iMethod).normalizationMedianSeconds = median(values(:,4));
        operations(iMethod).pipelineMedianSeconds = median(values(:,3)+values(:,4));
        operations(iMethod).internalMedianSeconds = median(values(:,6));
        operations(iMethod).boundaryResidualMedianSeconds = median(operations(iMethod).boundaryResidualSamplesSeconds);
        operations(iMethod).detectedCopyCountSamples = values(:,9).';
        operations(iMethod).detectedCopiedBytesSamples = values(:,10).';
    end
    [absoluteError,relativeError] = numericalError(outputs{iMethod},setup.reference);
    operations(iMethod).maximumAbsoluteError = absoluteError;
    operations(iMethod).maximumRelativeError = relativeError;
    operations(iMethod).correctnessPassed = relativeError <= options.errorTolerance;
end

existing = operations(ismember([operations.id],["dense-matrix","fft-extension"]));
fftw = operations([operations.id] == "fftw-allocating");
[existingTime,existingIndex] = min([existing.medianSeconds]);
validPrimary = operations(1:3);
validPrimary = validPrimary([validPrimary.correctnessPassed]);
[~,winnerIndex] = min([validPrimary.medianSeconds]);
workload = emptyWorkload;
workload.Nz = n;
workload.batchCount = batchCount;
workload.batchLabel = batchLabels(batchCount);
workload.dataType = dataType;
workload.transformType = transformType;
workload.direction = direction;
workload.nWarmups = options.nWarmups;
workload.nSamples = nSamples;
workload.operations = operations;
workload.fastestExistingMethod = existing(existingIndex).id;
workload.fastestExistingMedianSeconds = existingTime;
workload.winner = validPrimary(winnerIndex).id;
workload.fftwSpeedup = existingTime/fftw.medianSeconds;
workload.maximumRelativeError = max([operations.maximumRelativeError]);
workload.correctnessPassed = all([operations.correctnessPassed]);
workload.fftwEligible = fftw.correctnessPassed && workload.fftwSpeedup >= options.speedThreshold;
if ~workload.correctnessPassed
    error('FFTWR2RBenchmark:CorrectnessFailure','Correctness failed for Nz=%d, batch=%d, %s %s %s.',n,batchCount,dataType,transformType,direction);
end
end

function metrics = transformMetrics(transform)
% Access the plan only through a private test/benchmark bridge.
metrics = transform.backendMetrics;
end

function output = extensionTransform(input,coordinate,frequency,transformType,direction)
if transformType == "cosine" && direction == "forward"
    output = CosineTransformForward(coordinate,input,1);
elseif transformType == "cosine"
    output = CosineTransformBack(frequency,input,1);
elseif direction == "forward"
    output = SineTransformForward(coordinate,input,1,'both',0);
else
    output = SineTransformBack(frequency,input,1);
end
end

function eligibility = deriveEligibility(workloads,options)
eligibility = repmat(emptyEligibility,0,1);
for n = options.Nz
    for dataType = options.dataTypes
        for transformType = options.transformTypes
            for direction = ["forward","inverse"]
                selected = workloads([workloads.Nz] == n & [workloads.dataType] == dataType & [workloads.transformType] == transformType & [workloads.direction] == direction);
                [~,order] = sort([selected.batchCount]);
                selected = selected(order);
                passing = [selected.fftwEligible];
                intervals = repmat(struct('minimumBatchCount',0,'maximumBatchCount',0,'anchors',[]),0,1);
                startIndex = find(passing,1);
                while ~isempty(startIndex)
                    endIndex = startIndex;
                    while endIndex < numel(passing) && passing(endIndex+1)
                        endIndex = endIndex+1;
                    end
                    interval.minimumBatchCount = selected(startIndex).batchCount;
                    interval.maximumBatchCount = selected(endIndex).batchCount;
                    interval.anchors = [selected(startIndex:endIndex).batchCount];
                    intervals(end+1,1) = interval; %#ok<AGROW>
                    next = find(passing((endIndex+1):end),1);
                    if isempty(next)
                        startIndex = [];
                    else
                        startIndex = endIndex+next;
                    end
                end
                record = emptyEligibility;
                record.Nz = n;
                record.dataType = dataType;
                record.transformType = transformType;
                record.direction = direction;
                record.eligible = ~isempty(intervals);
                record.intervals = intervals;
                record.testedBatchCounts = [selected.batchCount];
                record.passingAnchors = [selected(passing).batchCount];
                eligibility(end+1,1) = record; %#ok<AGROW>
            end
        end
    end
end
end

function [forwardMatrix,backMatrix] = referenceMatrices(n,transformType)
if transformType == "cosine"
    forwardMatrix = CosineTransformForwardMatrix(n);
    backMatrix = CosineTransformBackMatrix(n);
else
    forwardMatrix = SineTransformForwardMatrix(n);
    backMatrix = SineTransformBackMatrix(n);
end
end

function output = zerosLike(sz,dataType)
output = zeros(sz);
if dataType == "complex"
    output = complex(output);
end
end

function [absoluteError,relativeError] = numericalError(actual,reference)
absoluteError = max(abs(actual-reference),[],'all');
relativeError = absoluteError/max(max(abs(reference),[],'all'),eps);
end

function labels = batchLabels(counts)
labels = strings(size(counts));
for iCount = 1:numel(counts)
    switch counts(iCount)
        case 1
            labels(iCount) = "scalar";
        case 8320
            labels(iCount) = "half-x 128x128";
        case 33024
            labels(iCount) = "half-x 256x256";
        case 131584
            labels(iCount) = "half-x 512x512";
        otherwise
            labels(iCount) = "batch " + counts(iCount);
    end
end
if isscalar(labels), labels = labels(1); end
end

function runDirectory = prepareRunDirectory(outputDirectory,runId,shouldWrite)
runDirectory = "";
if ~shouldWrite, return; end
runDirectory = fullfile(outputDirectory,runId);
if isfolder(runDirectory) || isfile(runDirectory)
    error('FFTWR2RBenchmark:OutputExists','Benchmark output already exists at %s.',runDirectory);
end
[created,message] = mkdir(runDirectory);
if ~created
    error('FFTWR2RBenchmark:OutputCreationFailed','Unable to create %s: %s',runDirectory,message);
end
end

function environment = collectEnvironment()
environment.matlabVersion = string(version);
environment.matlabRelease = string(version('-release'));
environment.architecture = string(computer('arch'));
environment.mexExtension = string(mexext);
environment.operatingSystem = string(system_dependent('getos'));
environment.processor = "unknown";
environment.machineModel = "unknown";
environment.physicalMemory = "unknown";
[status,hardwareText] = system('system_profiler SPHardwareDataType -json');
if status == 0
    hardware = jsondecode(hardwareText).SPHardwareDataType(1);
    environment.processor = optionalStructString(hardware,'chip_type',environment.processor);
    environment.machineModel = optionalStructString(hardware,'machine_model',environment.machineModel);
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

function value = optionalStructString(input,fieldName,defaultValue)
if isfield(input,fieldName), value = string(input.(fieldName)); else, value = defaultValue; end
end

function sources = collectSourceProvenance(sourceDirectory)
names = ["runFFTWR2RBenchmark.m","RealToRealTransform.m","fftw_r2r.cpp","fftw_backend_support.hpp"];
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
bytes = typecast(digest.digest(),'uint8');
hash = string(lower(reshape(dec2hex(bytes,2).',1,[])));
end

function writeResultArtifacts(result,runDirectory)
writeTextFile(fullfile(runDirectory,"real-to-real-benchmark.json"),string(jsonencode(replaceMissingStrings(result))));
writeTextFile(fullfile(runDirectory,"summary.md"),createMarkdownSummary(result));
end

function summary = createMarkdownSummary(result)
lines = strings(0,1);
lines(end+1) = "# FFTW real-to-real benchmark";
lines(end+1) = "";
lines(end+1) = "Status: **" + upper(result.status) + "**";
lines(end+1) = "";
lines(end+1) = "## Environment";
lines(end+1) = "";
lines(end+1) = "| Field | Value |";
lines(end+1) = "|---|---|";
lines(end+1) = markdownRow("MATLAB",result.environment.matlabVersion);
lines(end+1) = markdownRow("Architecture",result.environment.architecture);
lines(end+1) = markdownRow("Processor",result.environment.processor);
lines(end+1) = markdownRow("Memory",result.environment.physicalMemory);
lines(end+1) = markdownRow("FFTW",result.engine.version);
lines(end+1) = markdownRow("Library",result.engine.library);
if ~isempty(result.workloads)
    lines(end+1) = "";
    lines(end+1) = "## Complete-call winners";
    lines(end+1) = "";
    lines(end+1) = "| Nz | Batch | Type | Transform | Direction | Dense (ms) | Extension (ms) | FFTW (ms) | Speedup | Winner | Eligible | Max rel. error |";
    lines(end+1) = "|---:|---:|---|---|---|---:|---:|---:|---:|---|---|---:|";
    for workload = result.workloads'
        dense = workload.operations([workload.operations.id] == "dense-matrix");
        extension = workload.operations([workload.operations.id] == "fft-extension");
        fftw = workload.operations([workload.operations.id] == "fftw-allocating");
        lines(end+1) = sprintf("| %d | %d | %s | %s | %s | %.4f | %.4f | %.4f | %.3fx | %s | %s | %.3g |",workload.Nz,workload.batchCount,workload.dataType,workload.transformType,workload.direction,1e3*dense.medianSeconds,1e3*extension.medianSeconds,1e3*fftw.medianSeconds,workload.fftwSpeedup,workload.winner,yesNo(workload.fftwEligible),workload.maximumRelativeError); %#ok<AGROW>
    end
end
if ~isempty(result.eligibility)
    lines(end+1) = "";
    lines(end+1) = "## Bounded eligibility";
    lines(end+1) = "";
    lines(end+1) = "| Nz | Type | Transform | Direction | Eligible batch intervals |";
    lines(end+1) = "|---:|---|---|---|---|";
    for record = result.eligibility'
        intervalText = "none";
        if record.eligible
            values = arrayfun(@(interval) string(interval.minimumBatchCount) + "-" + string(interval.maximumBatchCount),record.intervals);
            intervalText = strjoin(values,", ");
        end
        lines(end+1) = sprintf("| %d | %s | %s | %s | %s |",record.Nz,record.dataType,record.transformType,record.direction,intervalText); %#ok<AGROW>
    end
end
lines(end+1) = "";
lines(end+1) = "## Timing boundaries";
lines(end+1) = "";
lines(end+1) = "- Complete-call time includes the MATLAB method or reference transform call.";
lines(end+1) = "- FFTW kernel time contains fftw_execute_r2r only.";
lines(end+1) = "- FFTW pipeline time contains kernel execution plus GL normalization and endpoint handling.";
lines(end+1) = "- Input generation, matrix construction, planning, preallocation, and diagnostic retrieval are excluded.";
if result.status == "failed" && ~isempty(result.failure)
    lines(end+1) = "";
    lines(end+1) = "## Failure";
    lines(end+1) = "";
    lines(end+1) = "- Stage: " + result.failure.stage;
    lines(end+1) = "- Identifier: " + result.failure.identifier;
    lines(end+1) = "- Message: " + replace(result.failure.message,"|","\|");
end
summary = strjoin(lines,newline);
end

function writeTextFile(path,text)
[fileId,message] = fopen(path,'w');
if fileId < 0
    error('FFTWR2RBenchmark:ArtifactOpenFailed','Unable to open %s: %s',path,message);
end
cleanup = onCleanup(@() fclose(fileId));
fprintf(fileId,'%s\n',text);
clear cleanup
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
row = "| " + replace(string(field),"|","\|") + " | " + replace(string(value),"|","\|") + " |";
end

function value = yesNo(condition)
if condition, value = "yes"; else, value = "no"; end
end

function operation = emptyOperation()
operation = struct('id',"",'totalSamplesSeconds',[],'medianSeconds',NaN,'allocationSamplesSeconds',[],'allocationMedianSeconds',NaN,'wrapSamplesSeconds',[],'wrapMedianSeconds',NaN,'kernelSamplesSeconds',[],'kernelMedianSeconds',NaN,'normalizationSamplesSeconds',[],'normalizationMedianSeconds',NaN,'pipelineSamplesSeconds',[],'pipelineMedianSeconds',NaN,'detachSamplesSeconds',[],'internalSamplesSeconds',[],'internalMedianSeconds',NaN,'boundaryResidualSamplesSeconds',[],'boundaryResidualMedianSeconds',NaN,'detectedCopyCountSamples',[],'detectedCopiedBytesSamples',[],'maximumAbsoluteError',NaN,'maximumRelativeError',NaN,'correctnessPassed',false);
end

function workload = emptyWorkload()
workload = struct('Nz',0,'batchCount',0,'batchLabel',"",'dataType',"",'transformType',"",'direction',"",'nWarmups',0,'nSamples',0,'operations',repmat(emptyOperation,0,1),'fastestExistingMethod',"",'fastestExistingMedianSeconds',NaN,'winner',"",'fftwSpeedup',NaN,'maximumRelativeError',NaN,'correctnessPassed',false,'fftwEligible',false);
end

function eligibility = emptyEligibility()
eligibility = struct('Nz',0,'dataType',"",'transformType',"",'direction',"",'eligible',false,'intervals',repmat(struct('minimumBatchCount',0,'maximumBatchCount',0,'anchors',[]),0,1),'testedBatchCounts',[],'passingAnchors',[]);
end
