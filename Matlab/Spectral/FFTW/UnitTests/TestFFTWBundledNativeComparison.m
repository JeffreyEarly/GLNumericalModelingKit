classdef TestFFTWBundledNativeComparison < matlab.unittest.TestCase
    methods (Test)
        function testMatchedComparisonReplayAndArtifacts(testCase)
            import matlab.unittest.fixtures.TemporaryFolderFixture

            fixture = testCase.applyFixture(TemporaryFolderFixture);
            [fftwDirectory,repositoryRoot] = TestFFTWBundledNativeComparison.repositoryPaths;
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
            issue38Path = fullfile(fixture.Folder,"issue38-smoke.json");
            TestFFTWBundledNativeComparison.writeIssue38Fixture(issue38Path);

            shouldBuild = exist('fftw_ownership_benchmark_bundled','file') ~= 3 || exist('fftw_ownership_benchmark_native','file') ~= 3 || exist('fftw_engine_benchmark_bundled','file') ~= 3 || exist('fftw_engine_benchmark_native','file') ~= 3;
            result = runFFTWBundledNativeComparison(outputDirectory=string(fixture.Folder),runId="success",issue38ArtifactPath=string(issue38Path),sizes=[16 8 4],gateSizes=[16 8 4],planners="estimate",alignmentModes=["matched","unaligned"],threadCount=1,plannerTimeLimitSeconds=1,nWarmups=1,nSamples=2,nSamplesLargest=2,shouldReplayIssue38=true,shouldBuild=shouldBuild);

            testCase.verifyEqual(result.status,"passed");
            testCase.verifyNumElements(result.engines,2);
            testCase.verifyNotEqual(result.engines(1).library,result.engines(2).library);
            testCase.verifyFalse(startsWith(result.engines(2).library,string(matlabroot)));
            testCase.verifyTrue(all([result.engines.alignmentMatchedAccepted]));
            testCase.verifyTrue(all([result.engines.alignmentMismatchRejected]));
            testCase.verifyTrue(all([result.engines.unalignedAccepted]));

            testCase.verifyNumElements(result.comparisons,2);
            testCase.verifyEqual([result.comparisons.status],["passed","passed"]);
            testCase.verifyEqual([result.comparisons.nWarmups],[1 1]);
            testCase.verifyEqual([result.comparisons.nSamples],[2 2]);
            configurations = [result.comparisons.configuration];
            testCase.verifyEqual([configurations.strategy],["guru-rank2","guru-rank2"]);
            testCase.verifyEqual(vertcat(configurations.transformOrder),[2 1; 2 1]);
            testCase.verifyEqual(unique([configurations.alignmentMode]),["matched","unaligned"]);
            for comparison = result.comparisons'
                testCase.verifyTrue(comparison.bundled.plan.wisdomClearedBeforePlanning);
                testCase.verifyTrue(comparison.native.plan.wisdomClearedBeforePlanning);
                testCase.verifyNumElements(comparison.matlab.forward.samplesSeconds,2);
                testCase.verifyNumElements(comparison.bundled.forward.totalMexSamplesSeconds,2);
                testCase.verifyNumElements(comparison.native.inverse.totalMexSamplesSeconds,2);
                testCase.verifyTrue(all(isfinite(comparison.bundled.forward.kernelSamplesSeconds)));
                testCase.verifyTrue(all(isfinite(comparison.native.forward.internalPipelineSamplesSeconds)));
                testCase.verifyGreaterThanOrEqual(comparison.bundled.forward.totalMexMedianSeconds,comparison.bundled.forward.kernelMedianSeconds);
                testCase.verifyGreaterThanOrEqual(comparison.native.forward.totalMexMedianSeconds,comparison.native.forward.kernelMedianSeconds);
                testCase.verifyTrue(comparison.bundled.correctnessPassed);
                testCase.verifyTrue(comparison.native.correctnessPassed);
                testCase.verifyTrue(comparison.bundled.zeroCopyPassed);
                testCase.verifyTrue(comparison.native.zeroCopyPassed);
                testCase.verifyTrue(comparison.bundled.destructiveInputChanged);
                testCase.verifyTrue(comparison.native.destructiveInputChanged);
                testCase.verifyEqual(comparison.bundledRelativeToNative.forwardTotalPercent,100*(comparison.bundled.forward.totalMexMedianSeconds/comparison.native.forward.totalMexMedianSeconds-1),'AbsTol',1e-12);
            end

            realElements = prod([16 8 4]);
            halfElements = (floor(16/2)+1)*8*4;
            testCase.verifyEqual(result.workloads.storage.realInputBytes,8*realElements);
            testCase.verifyEqual(result.workloads.storage.matlabFullSpectrumBytes,16*realElements);
            testCase.verifyEqual(result.workloads.storage.halfSpectrumBytes,16*halfElements);
            testCase.verifyLessThanOrEqual(result.workloads.bestBundled.maximumRelativeError,1e-12);
            testCase.verifyTrue(result.workloads.thresholds.correctnessPassed);
            testCase.verifyTrue(result.workloads.thresholds.storagePassed);
            testCase.verifyTrue(result.workloads.thresholds.zeroCopyPassed);

            testCase.verifyNumElements(result.historicalReplay,2);
            testCase.verifyEqual([result.historicalReplay.role],["best-bundled-screening","selected-native-winner"]);
            testCase.verifyTrue(all([result.historicalReplay.diagnosticOnly]));
            testCase.verifyTrue(all(arrayfun(@(replay) replay.nSamples == 2,result.historicalReplay)));
            testCase.verifyTrue(all(arrayfun(@(replay) replay.plan.wisdomClearedBeforePlanning,result.historicalReplay)));
            testCase.verifyTrue(all([result.historicalReplay.correctnessPassed]));

            jsonPath = fullfile(fixture.Folder,"success","bundled-native-comparison.json");
            markdownPath = fullfile(fixture.Folder,"success","summary.md");
            testCase.verifyTrue(isfile(jsonPath));
            testCase.verifyTrue(isfile(markdownPath));
            decoded = jsondecode(fileread(jsonPath));
            testCase.verifyEqual(string(decoded.status),"passed");
            testCase.verifyNumElements(decoded.comparisons,2);
            markdown = string(fileread(markdownPath));
            testCase.verifyTrue(contains(markdown,"## Environment"));
            testCase.verifyTrue(contains(markdown,"## Issue #38 discrepancy"));
            testCase.verifyTrue(contains(markdown,"## Matched configuration comparison"));
            testCase.verifyTrue(contains(markdown,"## Readiness by workload"));
            testCase.verifyTrue(contains(markdown,"## Criterion status"));
            testCase.verifyTrue(contains(markdown,"## Timing boundaries"));
            testCase.verifyEqual(string(fftw('planner')),previousPlanner);
            testCase.verifyEqual(fftw('dwisdom'),previousWisdom);
            testCase.verifyEqual(maxNumCompThreads,previousThreads);
            testCase.verifyEqual(rng,previousRandomState);
            clear pathCleanup directoryCleanup
        end

        function testFailureArtifactAndStateRestoration(testCase)
            import matlab.unittest.fixtures.TemporaryFolderFixture

            fixture = testCase.applyFixture(TemporaryFolderFixture);
            [fftwDirectory,repositoryRoot] = TestFFTWBundledNativeComparison.repositoryPaths;
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
            benchmark = @() runFFTWBundledNativeComparison(outputDirectory=string(fixture.Folder),runId="failure",sizes=[16 8 4],gateSizes=[16 8 4],planners="estimate",alignmentModes="unaligned",threadCount=1,plannerTimeLimitSeconds=1,nWarmups=0,nSamples=1,nSamplesLargest=1,errorTolerance=0,shouldReplayIssue38=false,shouldBuild=false);
            testCase.verifyError(benchmark,'FFTWBundledNative:NoValidEngineCandidate');

            jsonPath = fullfile(fixture.Folder,"failure","bundled-native-comparison.json");
            markdownPath = fullfile(fixture.Folder,"failure","summary.md");
            testCase.verifyTrue(isfile(jsonPath));
            testCase.verifyTrue(isfile(markdownPath));
            decoded = jsondecode(fileread(jsonPath));
            testCase.verifyEqual(string(decoded.status),"failed");
            testCase.verifyEqual(string(decoded.failure.identifier),"FFTWBundledNative:NoValidEngineCandidate");
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

        function writeIssue38Fixture(path)
            screening.forward.totalMexMedianSeconds = 1e-3;
            screening.forward.rawPipelineMedianSeconds = 8e-4;
            bundled.engine = "bundled-fftw";
            bundled.status = "passed";
            bundled.size = [16 8 4];
            bundled.transformOrder = [2 1];
            bundled.layout = "half-x";
            bundled.strategy = "guru-rank2";
            bundled.planner = "estimate";
            bundled.threads = 1;
            bundled.alignmentMode = "unaligned";
            bundled.screening = screening;
            native = bundled;
            native.engine = "native-fftw";
            native.screening.forward.totalMexMedianSeconds = 9e-4;
            native.screening.forward.rawPipelineMedianSeconds = 7e-4;
            artifact.candidates = [bundled native];
            artifact.workloads.size = [16 8 4];
            artifact.workloads.winner = native;
            [fileId,message] = fopen(path,'w');
            if fileId < 0, error('FFTWBundledNativeTest:FixtureOpenFailed','Unable to open fixture: %s',message); end
            cleanup = onCleanup(@() fclose(fileId));
            fprintf(fileId,'%s\n',jsonencode(artifact));
            clear cleanup
        end
    end
end
