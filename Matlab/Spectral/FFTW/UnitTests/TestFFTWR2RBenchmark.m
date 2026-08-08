classdef TestFFTWR2RBenchmark < matlab.unittest.TestCase
    methods (Test)
        function testSmokeArtifactsAndEligibility(testCase)
            import matlab.unittest.fixtures.TemporaryFolderFixture

            fixture = testCase.applyFixture(TemporaryFolderFixture);
            [fftwDirectory,repositoryRoot] = TestFFTWR2RBenchmark.repositoryPaths;
            previousDirectory = pwd;
            directoryCleanup = onCleanup(@() cd(previousDirectory));
            cd(repositoryRoot);
            previousPath = path;
            pathCleanup = onCleanup(@() path(previousPath));
            addpath(fullfile(repositoryRoot,'Matlab','Spectral'));
            addpath(fftwDirectory);
            if exist('fftw_r2r','file') ~= 3
                RealToRealTransform.makeMexFiles;
            end

            previousPlanner = string(fftw('planner'));
            previousWisdom = fftw('dwisdom');
            previousThreads = maxNumCompThreads;
            previousRandomState = rng;
            result = runFFTWR2RBenchmark(outputDirectory=string(fixture.Folder),runId="success",Nz=9,batchCounts=[1 4],dataTypes="real",transformTypes=["cosine","sine"],planner="estimate",threadCount=1,nWarmups=1,nSamples=2,nSamplesLargest=2,shouldBuild=false,requireCanonicalPlatform=false);

            testCase.verifyEqual(result.status,"passed");
            testCase.verifyNumElements(result.workloads,8);
            testCase.verifyNumElements(result.eligibility,4);
            for workload = result.workloads'
                testCase.verifyNumElements(workload.operations,4);
                testCase.verifyEqual([workload.operations.id],["dense-matrix","fft-extension","fftw-allocating","fftw-preallocated"]);
                testCase.verifyTrue(all(arrayfun(@(operation) numel(operation.totalSamplesSeconds) == 2,workload.operations)));
                testCase.verifyTrue(all(isfinite([workload.operations.medianSeconds])));
                testCase.verifyLessThanOrEqual(workload.maximumRelativeError,1e-12);
                fftwOperation = workload.operations([workload.operations.id] == "fftw-allocating");
                testCase.verifyNumElements(fftwOperation.kernelSamplesSeconds,2);
                testCase.verifyNumElements(fftwOperation.normalizationSamplesSeconds,2);
                testCase.verifyTrue(all(isfinite(fftwOperation.pipelineSamplesSeconds)));
                testCase.verifyEqual(workload.fftwEligible,workload.fftwSpeedup >= 1.10);
            end
            for record = result.eligibility'
                testCase.verifyEqual(record.testedBatchCounts,[1 4]);
                if record.eligible
                    testCase.verifyTrue(all(arrayfun(@(interval) interval.minimumBatchCount >= 1 && interval.maximumBatchCount <= 4,record.intervals)));
                end
            end

            jsonPath = fullfile(fixture.Folder,"success","real-to-real-benchmark.json");
            markdownPath = fullfile(fixture.Folder,"success","summary.md");
            testCase.verifyTrue(isfile(jsonPath));
            testCase.verifyTrue(isfile(markdownPath));
            decoded = jsondecode(fileread(jsonPath));
            testCase.verifyEqual(string(decoded.status),"passed");
            testCase.verifyNumElements(decoded.workloads,8);
            markdown = string(fileread(markdownPath));
            testCase.verifyTrue(contains(markdown,"## Environment"));
            testCase.verifyTrue(contains(markdown,"## Complete-call winners"));
            testCase.verifyTrue(contains(markdown,"## Bounded eligibility"));
            testCase.verifyTrue(contains(markdown,"## Timing boundaries"));
            testCase.verifyEqual(string(fftw('planner')),previousPlanner);
            testCase.verifyEqual(TestFFTWR2RBenchmark.canonicalWisdom(fftw('dwisdom')),TestFFTWR2RBenchmark.canonicalWisdom(previousWisdom));
            testCase.verifyEqual(maxNumCompThreads,previousThreads);
            testCase.verifyEqual(rng,previousRandomState);
            clear pathCleanup directoryCleanup
        end

        function testFailureArtifactAndStateRestoration(testCase)
            import matlab.unittest.fixtures.TemporaryFolderFixture

            fixture = testCase.applyFixture(TemporaryFolderFixture);
            [fftwDirectory,repositoryRoot] = TestFFTWR2RBenchmark.repositoryPaths;
            previousDirectory = pwd;
            directoryCleanup = onCleanup(@() cd(previousDirectory));
            cd(repositoryRoot);
            previousPath = path;
            pathCleanup = onCleanup(@() path(previousPath));
            addpath(fullfile(repositoryRoot,'Matlab','Spectral'));
            addpath(fftwDirectory);
            previousPlanner = string(fftw('planner'));
            previousWisdom = fftw('dwisdom');
            previousThreads = maxNumCompThreads;
            previousRandomState = rng;

            benchmark = @() runFFTWR2RBenchmark(outputDirectory=string(fixture.Folder),runId="failure",Nz=9,batchCounts=1,dataTypes="real",transformTypes="cosine",planner="estimate",threadCount=1,nWarmups=0,nSamples=1,nSamplesLargest=1,errorTolerance=0,shouldBuild=false,requireCanonicalPlatform=false);
            testCase.verifyError(benchmark,'FFTWR2RBenchmark:CorrectnessFailure');
            jsonPath = fullfile(fixture.Folder,"failure","real-to-real-benchmark.json");
            markdownPath = fullfile(fixture.Folder,"failure","summary.md");
            testCase.verifyTrue(isfile(jsonPath));
            testCase.verifyTrue(isfile(markdownPath));
            decoded = jsondecode(fileread(jsonPath));
            testCase.verifyEqual(string(decoded.status),"failed");
            testCase.verifyEqual(string(decoded.failure.identifier),"FFTWR2RBenchmark:CorrectnessFailure");
            testCase.verifyTrue(contains(string(fileread(markdownPath)),"## Failure"));
            testCase.verifyEqual(string(fftw('planner')),previousPlanner);
            testCase.verifyEqual(TestFFTWR2RBenchmark.canonicalWisdom(fftw('dwisdom')),TestFFTWR2RBenchmark.canonicalWisdom(previousWisdom));
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

        function wisdom = canonicalWisdom(wisdom)
            % FFTW may serialize an unchanged set of wisdom records in a
            % different order after import. Sort the complete record lines
            % so the test checks preserved content rather than text order.
            lines = strip(splitlines(string(wisdom)));
            wisdom = sort(lines(strlength(lines) > 0));
        end
    end
end
