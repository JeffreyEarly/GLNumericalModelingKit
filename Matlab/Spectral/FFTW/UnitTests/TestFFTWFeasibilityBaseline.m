classdef TestFFTWFeasibilityBaseline < matlab.unittest.TestCase
    methods (Test)
        function testCompleteSmokeRun(testCase)
            import matlab.unittest.fixtures.TemporaryFolderFixture

            fixture = testCase.applyFixture(TemporaryFolderFixture);
            [fftwDirectory,repositoryRoot] = TestFFTWFeasibilityBaseline.repositoryPaths();
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

            result = runFFTWFeasibilityBaseline(outputDirectory=string(fixture.Folder),sizes=[16 8 4],planner="estimate",nThreads=2,nWarmups=2,nSamples=2,nSamplesLargest=2,runId="smoke-success");

            testCase.verifyEqual(result.status,"passed");
            testCase.verifyNumElements(result.workloads,1);
            workload = result.workloads(1);
            testCase.verifyNumElements(workload.operations,7);
            testCase.verifyEqual(string({workload.operations.id}),["matlab_forward","matlab_inverse","mex_r2c_allocating","mex_r2c_preallocated","mex_c2r_allocating_preserving","mex_c2r_preallocated_preserving","mex_c2r_preallocated_destructive"]);
            testCase.verifyTrue(all(arrayfun(@(operation) numel(operation.sampleTimesSeconds) == 2,workload.operations)));
            testCase.verifyTrue(all(isfinite([workload.operations.medianSeconds])));
            testCase.verifyTrue(all([workload.operations.medianSeconds] > 0));
            testCase.verifyLessThanOrEqual(workload.correctness.maximumRelativeError,1e-12);
            testCase.verifyTrue(workload.correctness.preservingInputUnchanged);
            testCase.verifyTrue(workload.correctness.destructiveInputChanged);
            testCase.verifyTrue(workload.correctness.passed);

            expectedRealBytes = prod([16 8 4])*8;
            expectedFullSpectrumBytes = prod([16 8 4])*16;
            expectedHalfSpectrumBytes = 16*(8/2+1)*4*16;
            testCase.verifyEqual(workload.storage.realInputBytes,expectedRealBytes);
            testCase.verifyEqual(workload.storage.matlabFullSpectrumBytes,expectedFullSpectrumBytes);
            testCase.verifyEqual(workload.storage.mexHalfSpectrumBytes,expectedHalfSpectrumBytes);
            testCase.verifyEqual(workload.storage.halfToFullSpectrumStorageRatio,expectedHalfSpectrumBytes/expectedFullSpectrumBytes);

            jsonPath = fullfile(result.artifacts.directory,result.artifacts.json);
            markdownPath = fullfile(result.artifacts.directory,result.artifacts.markdown);
            testCase.verifyTrue(isfile(jsonPath));
            testCase.verifyTrue(isfile(markdownPath));
            jsonText = fileread(jsonPath);
            decoded = jsondecode(jsonText);
            testCase.verifyEqual(string(decoded.status),"passed");
            testCase.verifyTrue(contains(jsonText,'"workloads": ['));
            testCase.verifyNumElements(decoded.workloads,1);
            markdown = string(fileread(markdownPath));
            testCase.verifyTrue(contains(markdown,"## Median transform-call time"));
            testCase.verifyTrue(contains(markdown,"## Correctness"));
            testCase.verifyTrue(contains(markdown,"## Array storage"));
            testCase.verifyTrue(contains(markdown,"## Allocation model"));

            testCase.verifyEqual(string(fftw('planner')),previousPlanner);
            testCase.verifyEqual(fftw('dwisdom'),previousWisdom);
            testCase.verifyEqual(maxNumCompThreads,previousThreads);
            testCase.verifyEqual(rng,previousRandomState);

            clear pathCleanup directoryCleanup
        end

        function testFailureWritesArtifactAndRestoresState(testCase)
            import matlab.unittest.fixtures.TemporaryFolderFixture

            fixture = testCase.applyFixture(TemporaryFolderFixture);
            [fftwDirectory,repositoryRoot] = TestFFTWFeasibilityBaseline.repositoryPaths();
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

            benchmark = @() runFFTWFeasibilityBaseline(outputDirectory=string(fixture.Folder),sizes=[16 8 4],planner="estimate",nThreads=2,nWarmups=1,nSamples=1,nSamplesLargest=1,errorTolerance=0,runId="smoke-failure");
            testCase.verifyError(benchmark,'FFTWBenchmark:CorrectnessFailure');

            jsonPath = fullfile(fixture.Folder,'smoke-failure','benchmark.json');
            markdownPath = fullfile(fixture.Folder,'smoke-failure','summary.md');
            testCase.verifyTrue(isfile(jsonPath));
            testCase.verifyTrue(isfile(markdownPath));
            decoded = jsondecode(fileread(jsonPath));
            testCase.verifyEqual(string(decoded.status),"failed");
            testCase.verifyEqual(string(decoded.failure.identifier),"FFTWBenchmark:CorrectnessFailure");
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
