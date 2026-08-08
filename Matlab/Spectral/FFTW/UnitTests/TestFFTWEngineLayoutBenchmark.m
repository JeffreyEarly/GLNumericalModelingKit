classdef TestFFTWEngineLayoutBenchmark < matlab.unittest.TestCase
    methods (Test)
        function testCompleteSmokeRun(testCase)
            import matlab.unittest.fixtures.TemporaryFolderFixture

            fixture = testCase.applyFixture(TemporaryFolderFixture);
            [fftwDirectory,repositoryRoot] = TestFFTWEngineLayoutBenchmark.repositoryPaths();
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

            result = runFFTWEngineLayoutBenchmark(outputDirectory=string(fixture.Folder),runId="smoke-success",sizes=[16 8 4],gateSizes=[16 8 4],transformOrders=[1 2; 2 1; 3 1],planners="estimate",threadCounts=1,alignmentModes=["matched","unaligned"],nScreeningWarmups=1,nScreeningSamples=1,nWarmups=1,nSamples=2,nSamplesLargest=2,fallbackPolicy="always",shouldBuildNative=false);

            testCase.verifyEqual(result.status,"passed");
            testCase.verifyNumElements(result.workloads,1);
            testCase.verifyTrue(all(ismember(["bundled-fftw","accelerate-vdsp"],unique([result.candidates.engine]))));
            testCase.verifyTrue(all(ismember(["half-x","half-y"],unique([result.candidates.layout]))));
            testCase.verifyTrue(all(ismember(["guru-rank2","staged-r2c-c2c","vdsp-2d","vdsp-staged"],unique([result.candidates.strategy]))));
            testCase.verifyTrue(any(arrayfun(@(candidate) isequal(candidate.transformOrder,[3 1]),result.candidates)));
            testCase.verifyTrue(all(ismember(["matched","unaligned"],unique([result.candidates([result.candidates.engine] == "bundled-fftw").alignmentMode]))));

            passedCandidates = result.candidates([result.candidates.status] == "passed");
            testCase.verifyNotEmpty(passedCandidates);
            testCase.verifyTrue(all(isfinite(arrayfun(@(candidate) candidate.screening.forward.totalMexMedianSeconds,passedCandidates))));
            testCase.verifyTrue(all(isfinite(arrayfun(@(candidate) candidate.screening.forward.rawPipelineMedianSeconds,passedCandidates))));
            testCase.verifyTrue(all([passedCandidates.maximumRelativeError] <= 1e-12));
            testCase.verifyTrue(all([passedCandidates.destructiveInputChanged]));

            workload = result.workloads(1);
            testCase.verifyEqual(workload.finalMeasurement.nWarmups,1);
            testCase.verifyEqual(workload.finalMeasurement.nSamples,2);
            testCase.verifyNumElements(workload.finalMeasurement.matlabForwardSamplesSeconds,2);
            testCase.verifyNumElements(workload.finalMeasurement.forward.allocating.totalMexSamplesSeconds,2);
            testCase.verifyNumElements(workload.finalMeasurement.forward.preallocated.totalMexSamplesSeconds,2);
            testCase.verifyNumElements(workload.finalMeasurement.inverse.totalMexSamplesSeconds,2);
            testCase.verifyLessThanOrEqual(workload.finalMeasurement.maximumRelativeError,1e-12);
            testCase.verifyTrue(workload.finalMeasurement.destructiveInputChanged);
            testCase.verifyLessThanOrEqual(workload.finalMeasurement.forward.rawPipelineMedianSeconds,1.5*workload.finalMeasurement.forward.totalMexMedianSeconds);

            expectedRealBytes = prod([16 8 4])*8;
            compressedDimension = workload.winner.transformOrder(2);
            expectedComplexSize = [16 8 4];
            expectedComplexSize(compressedDimension) = expectedComplexSize(compressedDimension)/2+1;
            expectedHalfBytes = prod(expectedComplexSize)*16;
            testCase.verifyEqual(workload.storage.realInputBytes,expectedRealBytes);
            testCase.verifyEqual(workload.storage.matlabFullSpectrumBytes,2*expectedRealBytes);
            testCase.verifyEqual(workload.storage.halfSpectrumBytes,expectedHalfBytes);
            testCase.verifyEqual(workload.storage.halfToFullSpectrumStorageRatio,expectedHalfBytes/(2*expectedRealBytes));

            jsonPath = fullfile(result.artifacts.directory,result.artifacts.json);
            markdownPath = fullfile(result.artifacts.directory,result.artifacts.markdown);
            testCase.verifyTrue(isfile(jsonPath));
            testCase.verifyTrue(isfile(markdownPath));
            decoded = jsondecode(fileread(jsonPath));
            testCase.verifyEqual(string(decoded.status),"passed");
            testCase.verifyNumElements(decoded.workloads,1);
            testCase.verifyNumElements(decoded.candidates,numel(result.candidates));
            markdown = string(fileread(markdownPath));
            testCase.verifyTrue(contains(markdown,"## Environment"));
            testCase.verifyTrue(contains(markdown,"## Engines"));
            testCase.verifyTrue(contains(markdown,"## Fastest valid configuration by workload"));
            testCase.verifyTrue(contains(markdown,"## Criterion status"));
            testCase.verifyTrue(contains(markdown,"## Timing boundaries"));

            [matchedAccepted,mismatchRejected] = fftw_engine_benchmark_bundled('alignmentSelfTest');
            testCase.verifyTrue(matchedAccepted);
            testCase.verifyTrue(mismatchRejected);
            testCase.verifyEqual(string(fftw('planner')),previousPlanner);
            testCase.verifyEqual(fftw('dwisdom'),previousWisdom);
            testCase.verifyEqual(maxNumCompThreads,previousThreads);
            testCase.verifyEqual(rng,previousRandomState);
            clear pathCleanup directoryCleanup
        end

        function testFailureWritesPartialArtifactAndRestoresState(testCase)
            import matlab.unittest.fixtures.TemporaryFolderFixture

            fixture = testCase.applyFixture(TemporaryFolderFixture);
            [fftwDirectory,repositoryRoot] = TestFFTWEngineLayoutBenchmark.repositoryPaths();
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
            benchmark = @() runFFTWEngineLayoutBenchmark(outputDirectory=string(fixture.Folder),runId="smoke-failure",sizes=[16 8 4],gateSizes=[16 8 4],transformOrders=[1 2],fftwStrategies="guru-rank2",planners="estimate",threadCounts=1,alignmentModes="matched",nScreeningWarmups=0,nScreeningSamples=1,nWarmups=0,nSamples=1,nSamplesLargest=1,errorTolerance=0,fallbackPolicy="never",shouldBuild=false);
            testCase.verifyError(benchmark,'FFTWEngineBenchmark:NoValidCandidate');

            jsonPath = fullfile(fixture.Folder,"smoke-failure","engine-layout-benchmark.json");
            markdownPath = fullfile(fixture.Folder,"smoke-failure","summary.md");
            testCase.verifyTrue(isfile(jsonPath));
            testCase.verifyTrue(isfile(markdownPath));
            decoded = jsondecode(fileread(jsonPath));
            testCase.verifyEqual(string(decoded.status),"failed");
            testCase.verifyEqual(string(decoded.failure.identifier),"FFTWEngineBenchmark:NoValidCandidate");
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
