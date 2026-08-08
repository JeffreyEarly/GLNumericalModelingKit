classdef TestFFTWMexOwnershipBenchmark < matlab.unittest.TestCase
    methods (Test)
        function testSmokeBenchmarkAndArtifacts(testCase)
            import matlab.unittest.fixtures.TemporaryFolderFixture

            fixture = testCase.applyFixture(TemporaryFolderFixture);
            [fftwDirectory,repositoryRoot] = TestFFTWMexOwnershipBenchmark.repositoryPaths();
            previousDirectory = pwd;
            directoryCleanup = onCleanup(@() cd(previousDirectory));
            cd(repositoryRoot);
            previousPath = path;
            pathCleanup = onCleanup(@() path(previousPath));
            addpath(fftwDirectory);

            previousPlanner = string(fftw('planner'));
            previousWisdom = fftw('dwisdom');
            previousThreads = maxNumCompThreads;
            previousRandomState = rng;
            result = runFFTWMexOwnershipBenchmark(outputDirectory=string(fixture.Folder),runId="smoke",sizes=[16 8 4],gateSizes=[16 8 4],nWarmups=1,nSamples=2,nSamplesLargest=2,nPlanMemorySamples=1,shouldUseIssue38Winners=false,shouldBuild=true,shouldBuildNative=false);

            testCase.verifyEqual(result.status,"passed");
            testCase.verifyNumElements(result.workloads,1);
            workload = result.workloads(1);
            testCase.verifyEqual(workload.nWarmups,1);
            testCase.verifyEqual(workload.nSamples,2);
            testCase.verifyNumElements(workload.operations,11);
            expectedIds = ["forward-factory-array","forward-caller-direct","forward-caller-complex","forward-matlab-buffer","forward-fftw-buffer","inverse-preserving-allocating","inverse-preserving-preallocated","inverse-destructive-unique","inverse-destructive-aliased","inverse-destructive-complex","chain-fftw-forward-destructive-inverse"];
            testCase.verifyEqual([workload.operations.id],expectedIds);
            testCase.verifyTrue(all(arrayfun(@(operation) numel(operation.totalMexSamplesSeconds) == 2,workload.operations)));
            testCase.verifyTrue(all(isfinite([workload.operations.medianSeconds])));
            testCase.verifyLessThanOrEqual(workload.correctness.maximumRelativeError,1e-12);
            testCase.verifyTrue(workload.semantics.preservingInputUnchanged);
            testCase.verifyTrue(workload.semantics.destructiveInputChanged);
            testCase.verifyTrue(workload.semantics.aliasPreserved);

            operations = workload.operations;
            matlabBuffer = operations([operations.id] == "forward-matlab-buffer");
            fftwBuffer = operations([operations.id] == "forward-fftw-buffer");
            preserving = operations([operations.id] == "inverse-preserving-preallocated");
            destructive = operations([operations.id] == "inverse-destructive-unique");
            aliased = operations([operations.id] == "inverse-destructive-aliased");
            chain = operations([operations.id] == "chain-fftw-forward-destructive-inverse");
            testCase.verifyTrue(matlabBuffer.returnedPointerPreserved);
            testCase.verifyTrue(fftwBuffer.returnedPointerPreserved);
            testCase.verifyTrue(all(fftwBuffer.detectedCopiedBytesSamples == 0));
            testCase.verifyTrue(all(preserving.explicitCopyCountSamples == 1));
            testCase.verifyTrue(all(preserving.explicitCopiedBytesSamples == workload.storage.halfSpectrumBytes));
            testCase.verifyTrue(preserving.returnedPointerPreserved);
            testCase.verifyEqual(preserving.callerPreallocatedBytes,workload.storage.realArrayBytes);
            testCase.verifyFalse(any(destructive.inputDetached));
            testCase.verifyTrue(all(destructive.explicitCopiedBytesSamples == 0));
            testCase.verifyEqual(destructive.callerPreallocatedBytes,workload.storage.halfSpectrumBytes+workload.storage.realArrayBytes);
            testCase.verifyTrue(all(aliased.inputDetached));
            testCase.verifyTrue(chain.handoffPointerPreserved);
            testCase.verifyTrue(workload.lifetime.balanced);
            testCase.verifyEqual(workload.lifetime.created,6);
            testCase.verifyEqual(workload.lifetime.freed,6);
            testCase.verifyEqual(workload.lifetime.outstanding,0);
            testCase.verifyTrue(all([result.engines.alignmentMatchedAccepted]));
            testCase.verifyTrue(all([result.engines.alignmentMismatchRejected]));
            testCase.verifyTrue(all([result.engines.unalignedAccepted]));

            realElements = prod([16 8 4]);
            halfXElements = (16/2+1)*8*4;
            halfYElements = 16*(8/2+1)*4;
            testCase.verifyEqual(workload.storage.realArrayBytes,8*realElements);
            testCase.verifyEqual(workload.storage.matlabFullSpectrumBytes,16*realElements);
            testCase.verifyEqual(workload.storage.halfSpectrumBytes,16*halfXElements);
            testCase.verifyEqual(workload.storage.standaloneLegacy.totalKnownBytes,32*halfXElements);
            expectedWaveVortexFFTWBytes = 24*realElements + 48*(halfXElements+2*halfYElements);
            testCase.verifyEqual(workload.storage.waveVortexFFTW.totalKnownBytes,expectedWaveVortexFFTWBytes);
            testCase.verifyTrue(workload.plan.memoryEstimate.isEstimate);

            jsonPath = fullfile(fixture.Folder,"smoke","ownership-benchmark.json");
            markdownPath = fullfile(fixture.Folder,"smoke","summary.md");
            testCase.verifyTrue(isfile(jsonPath));
            testCase.verifyTrue(isfile(markdownPath));
            decoded = jsondecode(fileread(jsonPath));
            testCase.verifyEqual(string(decoded.status),"passed");
            testCase.verifyNumElements(decoded.workloads.operations,11);
            markdown = string(fileread(markdownPath));
            testCase.verifyTrue(contains(markdown,"## Ownership timing"));
            testCase.verifyTrue(contains(markdown,"## Zero-copy and lifetime evidence"));
            testCase.verifyTrue(contains(markdown,"## Persistent storage"));
            testCase.verifyTrue(contains(markdown,"## Recommendation"));
            testCase.verifyEqual(string(fftw('planner')),previousPlanner);
            testCase.verifyEqual(fftw('dwisdom'),previousWisdom);
            testCase.verifyEqual(maxNumCompThreads,previousThreads);
            testCase.verifyEqual(rng,previousRandomState);
            clear pathCleanup directoryCleanup
        end

        function testFailureArtifactAndStateRestoration(testCase)
            import matlab.unittest.fixtures.TemporaryFolderFixture

            fixture = testCase.applyFixture(TemporaryFolderFixture);
            [fftwDirectory,repositoryRoot] = TestFFTWMexOwnershipBenchmark.repositoryPaths();
            previousDirectory = pwd;
            directoryCleanup = onCleanup(@() cd(previousDirectory));
            cd(repositoryRoot);
            previousPath = path;
            pathCleanup = onCleanup(@() path(previousPath));
            addpath(fftwDirectory);

            previousPlanner = string(fftw('planner'));
            previousWisdom = fftw('dwisdom');
            previousThreads = maxNumCompThreads;
            previousRandomState = rng;
            benchmark = @() runFFTWMexOwnershipBenchmark(outputDirectory=string(fixture.Folder),runId="failure",sizes=[16 8 4],gateSizes=[16 8 4],nWarmups=0,nSamples=1,nSamplesLargest=1,nPlanMemorySamples=1,errorTolerance=0,shouldUseIssue38Winners=false,shouldBuild=false,shouldBuildNative=false);
            testCase.verifyError(benchmark,'FFTWMexOwnership:WorkloadValidationFailed');

            jsonPath = fullfile(fixture.Folder,"failure","ownership-benchmark.json");
            markdownPath = fullfile(fixture.Folder,"failure","summary.md");
            testCase.verifyTrue(isfile(jsonPath));
            testCase.verifyTrue(isfile(markdownPath));
            decoded = jsondecode(fileread(jsonPath));
            testCase.verifyEqual(string(decoded.status),"failed");
            testCase.verifyEqual(string(decoded.failure.identifier),"FFTWMexOwnership:WorkloadValidationFailed");
            testCase.verifyTrue(contains(string(fileread(markdownPath)),"## Failure"));
            testCase.verifyEqual(string(fftw('planner')),previousPlanner);
            testCase.verifyEqual(fftw('dwisdom'),previousWisdom);
            testCase.verifyEqual(maxNumCompThreads,previousThreads);
            testCase.verifyEqual(rng,previousRandomState);
            clear pathCleanup directoryCleanup
        end
    end

    methods (Static, Access=private)
        function [fftwDirectory,repositoryRoot] = repositoryPaths()
            unitTestDirectory = fileparts(mfilename('fullpath'));
            fftwDirectory = fileparts(unitTestDirectory);
            repositoryRoot = fileparts(fileparts(fileparts(fftwDirectory)));
        end
    end
end
