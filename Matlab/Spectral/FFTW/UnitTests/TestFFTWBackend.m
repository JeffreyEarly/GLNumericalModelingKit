classdef TestFFTWBackend < matlab.unittest.TestCase
    methods (TestClassSetup)
        function configurePathsAndBackend(testCase)
            testDirectory = fileparts(mfilename('fullpath'));
            fftwDirectory = fileparts(testDirectory);
            spectralDirectory = fileparts(fftwDirectory);
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(testDirectory));
            addpath(spectralDirectory);
            addpath(fftwDirectory);
            capabilities = FFTWBackend.build();
            testCase.assertTrue(capabilities.isComplete,capabilities.build.reason.message);
        end
    end

    methods (Test)
        function testCapabilitySchemaAndSelfTests(testCase)
            capabilities = FFTWBackend.capabilities();
            testCase.verifyEqual(capabilities.schemaVersion,"1.0.0");
            testCase.verifyEqual(capabilities.status,"available");
            testCase.verifyTrue(capabilities.isAvailable);
            testCase.verifyTrue(capabilities.isComplete);
            testCase.verifyEqual(capabilities.platform.matlabRelease,"2026a");
            testCase.verifyEqual(capabilities.platform.architecture,"maca64");
            testCase.verifyTrue(capabilities.platform.isSupported);
            testCase.verifyEqual(capabilities.provider.id,"matlab-bundled");
            testCase.verifyEqual(capabilities.provider.type,"bundled");
            testCase.verifyTrue(capabilities.provider.isDefault);
            testCase.verifyTrue(capabilities.provider.identityValidated);
            testCase.verifyTrue(capabilities.compiler.isAvailable);
            testCase.verifyTrue(capabilities.library.exists);
            testCase.verifyTrue(capabilities.library.identityValidated);
            testCase.verifyEqual(capabilities.library.version,"fftw-3.3.8");
            testCase.verifyEqual(capabilities.library.resolvedPath,capabilities.library.expectedPath);
            testCase.verifyTrue(capabilities.modules.r2c.isAvailable);
            testCase.verifyTrue(capabilities.modules.r2r.isAvailable);
            testCase.verifyTrue(endsWith(capabilities.modules.r2c.path,"fftw_r2c."+mexext));
            testCase.verifyTrue(endsWith(capabilities.modules.r2r.path,"fftw_r2r."+mexext));
            for name = ["r2c","c2r","dct1","dst1"]
                feature = capabilities.features.(name);
                testCase.verifyTrue(feature.isAvailable);
                testCase.verifyTrue(feature.selfTestPassed);
                testCase.verifyTrue(feature.alignmentPassed);
                testCase.verifyLessThanOrEqual(feature.maximumRelativeError,1e-12);
                testCase.verifyEqual(feature.reason.code,"");
            end
            testCase.verifyFalse(capabilities.build.attempted);
            testCase.verifyFalse(capabilities.build.isRequired);
            testCase.verifyTrue(capabilities.build.isPossible);
            testCase.verifyNumElements(capabilities.eligibility.realToReal.records,40);
        end

        function testEligibilityMatchesCanonicalArtifacts(testCase)
            capabilities = FFTWBackend.capabilities();
            repositoryRoot = TestFFTWBackend.repositoryRoot;
            issue41Path = fullfile(repositoryRoot,capabilities.eligibility.horizontal.sourceArtifact);
            issue43Path = fullfile(repositoryRoot,capabilities.eligibility.realToReal.sourceArtifact);
            issue41 = jsondecode(fileread(issue41Path));
            issue43 = jsondecode(fileread(issue43Path));
            testCase.verifyTrue(issue41.readiness.ready);
            testCase.verifyEqual(capabilities.eligibility.horizontal.gateSizes,issue41.configuration.gateSizes);
            testCase.verifyEqual(capabilities.eligibility.horizontal.thresholds.rawForwardSpeedup,issue41.configuration.thresholds.rawForwardSpeedRatio);
            testCase.verifyEqual(capabilities.eligibility.horizontal.thresholds.completeForwardSpeedup,issue41.configuration.thresholds.totalForwardSpeedRatio);
            testCase.verifyEqual(capabilities.eligibility.horizontal.thresholds.destructiveInverseSpeedRatio,issue41.configuration.thresholds.inverseSpeedRatio);

            records = capabilities.eligibility.realToReal.records;
            testCase.verifyNumElements(issue43.eligibility,numel(records));
            for record = records'
                artifact = issue43.eligibility([issue43.eligibility.Nz] == record.Nz & string({issue43.eligibility.dataType}) == record.dataType & string({issue43.eligibility.transformType}) == record.transformType & string({issue43.eligibility.direction}) == record.direction);
                testCase.verifyNumElements(artifact,1);
                testCase.verifyEqual(record.eligible,artifact.eligible);
                testCase.verifyEqual(record.testedBatchCounts,artifact.testedBatchCounts(:).');
                if record.eligible
                    testCase.verifyEqual(record.intervals.minimumBatchCount,artifact.intervals.minimumBatchCount);
                    testCase.verifyEqual(record.intervals.maximumBatchCount,artifact.intervals.maximumBatchCount);
                else
                    testCase.verifyEmpty(record.intervals);
                end
            end

            complexDct65 = records([records.Nz] == 65 & [records.dataType] == "complex" & [records.transformType] == "cosine" & [records.direction] == "forward");
            testCase.verifyFalse(TestFFTWBackend.isEligible(complexDct65,8320));
            testCase.verifyTrue(TestFFTWBackend.isEligible(complexDct65,33024));
            testCase.verifyTrue(TestFFTWBackend.isEligible(complexDct65,100000));
            testCase.verifyFalse(TestFFTWBackend.isEligible(complexDct65,131585));
            testCase.verifyFalse(any([records.Nz] == 64));
        end

        function testStructuredUnavailableAndPartialResults(testCase)
            context = TestableFFTWBackend.context();
            context.release = "2025b";
            unavailable = TestableFFTWBackend.inspect(context);
            testCase.verifyEqual(unavailable.status,"unavailable");
            testCase.verifyEqual(unavailable.reason.code,"unsupported-release");
            testCase.verifyFalse(unavailable.build.isPossible);

            context = TestableFFTWBackend.context();
            context.architecture = "glnxa64";
            context.provider.supportedArchitectures = "maca64";
            unavailable = TestableFFTWBackend.inspect(context);
            testCase.verifyEqual(unavailable.reason.code,"unsupported-architecture");

            context = TestableFFTWBackend.context();
            context.failureInjection = "library-missing";
            unavailable = TestableFFTWBackend.inspect(context);
            testCase.verifyEqual(unavailable.reason.code,"bundled-library-missing");

            context = TestableFFTWBackend.context();
            context.failureInjection = "fftw_r2c-missing";
            partial = TestableFFTWBackend.inspect(context);
            testCase.verifyEqual(partial.status,"partial");
            testCase.verifyTrue(partial.isAvailable);
            testCase.verifyFalse(partial.isComplete);
            testCase.verifyFalse(partial.features.r2c.isAvailable);
            testCase.verifyTrue(partial.features.dct1.isAvailable);
            testCase.verifyEqual(partial.features.r2c.reason.code,"mex-missing");

            context = TestableFFTWBackend.context();
            context.failureInjection = "fftw_r2r-missing";
            partial = TestableFFTWBackend.inspect(context);
            testCase.verifyTrue(partial.features.r2c.isAvailable);
            testCase.verifyFalse(partial.features.dct1.isAvailable);
            testCase.verifyEqual(partial.features.dct1.reason.code,"mex-missing");

            expectations = {
                "unexpected-module-path","unexpected-module-path"
                "missing-symbol","missing-symbol"
                "library-mismatch","library-mismatch"
                "alignment-failed","alignment-failed"
                "r2c-self-test-failed","self-test-failed"
                "dst1-self-test-failed","self-test-failed"
                };
            for iExpectation = 1:size(expectations,1)
                context = TestableFFTWBackend.context();
                context.failureInjection = expectations{iExpectation,1};
                result = TestableFFTWBackend.inspect(context);
                codes = [result.reason.code result.features.r2c.reason.code result.features.c2r.reason.code result.features.dct1.reason.code result.features.dst1.reason.code result.modules.r2c.reason.code result.modules.r2r.reason.code];
                testCase.verifyTrue(any(codes == expectations{iExpectation,2}),"Missing structured failure for "+expectations{iExpectation,1});
            end
        end

        function testProviderNeutralOrchestrationAndStateRestoration(testCase)
            previousDirectory = pwd;
            previousPath = path;
            previousPlanner = string(fftw('planner'));
            previousWisdom = fftw('dwisdom');
            previousThreads = maxNumCompThreads;
            previousRandom = rng;

            context = TestableFFTWBackend.context();
            context.provider.id = "test-provider";
            context.provider.type = "test";
            context.provider.origin = "injected test provider";
            result = TestableFFTWBackend.inspect(context);
            testCase.verifyEqual(result.status,"available");
            testCase.verifyEqual(result.provider.id,"test-provider");
            testCase.verifyEqual(result.provider.type,"test");
            testCase.verifyTrue(result.provider.identityValidated);

            context.failureInjection = "alignment-failed";
            TestableFFTWBackend.inspect(context);
            testCase.verifyEqual(pwd,previousDirectory);
            testCase.verifyEqual(path,previousPath);
            testCase.verifyEqual(string(fftw('planner')),previousPlanner);
            testCase.verifyEqual(TestFFTWBackend.canonicalWisdom(fftw('dwisdom')),TestFFTWBackend.canonicalWisdom(previousWisdom));
            testCase.verifyEqual(maxNumCompThreads,previousThreads);
            testCase.verifyEqual(rng,previousRandom);
        end

        function testCapabilityQueryPreservesLiveTransform(testCase)
            input = reshape(sin((1:(8*6*3))/7),[8 6 3]);
            transform = RealToComplexTransform(size(input),dims=[2 1],planner="estimate",nCores=1);
            cleanup = onCleanup(@() delete(transform));
            expected = transform.transformForward(input);
            lifetimeBefore = fftw_r2c('lifetime');
            testCase.verifyTrue(mislocked('fftw_r2c'));
            capabilities = FFTWBackend.capabilities();
            actual = transform.transformForward(input);
            lifetimeAfter = fftw_r2c('lifetime');
            testCase.verifyEqual(capabilities.status,"available");
            testCase.verifyEqual(lifetimeAfter(3),lifetimeBefore(3));
            testCase.verifyEqual(actual,expected);
            blockedBuild = FFTWBackend.build();
            testCase.verifyFalse(blockedBuild.build.succeeded);
            testCase.verifyEqual(blockedBuild.build.reason.code,"module-in-use");
            clear cleanup
        end

        function testBuildFromAnyWorkingDirectory(testCase)
            import matlab.unittest.fixtures.TemporaryFolderFixture
            fixture = testCase.applyFixture(TemporaryFolderFixture);
            previousDirectory = pwd;
            directoryCleanup = onCleanup(@() cd(previousDirectory));

            cd(TestFFTWBackend.repositoryRoot);
            rootBuild = FFTWBackend.build();
            testCase.verifyTrue(rootBuild.build.attempted);
            testCase.verifyTrue(rootBuild.build.succeeded,rootBuild.build.reason.message);
            testCase.verifyTrue(rootBuild.build.installed);
            testCase.verifyTrue(rootBuild.isComplete);
            testCase.verifyNotEqual(rootBuild.build.startedAtUTC,"");
            testCase.verifyNotEqual(rootBuild.build.completedAtUTC,"");
            testCase.verifyGreaterThan(rootBuild.build.durationSeconds,0);

            cd(fixture.Folder);
            unrelatedBuild = FFTWBackend.build();
            testCase.verifyTrue(unrelatedBuild.build.succeeded,unrelatedBuild.build.reason.message);
            testCase.verifyTrue(unrelatedBuild.isComplete);
            clear directoryCleanup
        end

        function testLockedModuleAndTransactionalRollback(testCase)
            context = TestableFFTWBackend.context();
            modulePaths = TestFFTWBackend.modulePaths(context);
            baselineHashes = arrayfun(@TestFFTWBackend.fileHash,modulePaths);

            transform = RealToComplexTransform([8 6 3],dims=[2 1],planner="estimate",nCores=1);
            transformCleanup = onCleanup(@() delete(transform));
            lockedResult = TestableFFTWBackend.buildUsing(context);
            testCase.verifyFalse(lockedResult.build.succeeded);
            testCase.verifyTrue(contains(lockedResult.build.reason.identifier,"ModuleInUse"));
            testCase.verifyEqual(arrayfun(@TestFFTWBackend.fileHash,modulePaths),baselineHashes);
            clear transformCleanup

            context.failureInjection = "compile-failed";
            compileFailure = TestableFFTWBackend.buildUsing(context);
            testCase.verifyFalse(compileFailure.build.succeeded);
            testCase.verifyEqual(compileFailure.build.reason.code,"compile-failed");
            testCase.verifyTrue(contains(compileFailure.build.reason.identifier,"CompileFailure"));
            testCase.verifyEqual(arrayfun(@TestFFTWBackend.fileHash,modulePaths),baselineHashes);

            context.failureInjection = "install-failed";
            installFailure = TestableFFTWBackend.buildUsing(context);
            testCase.verifyFalse(installFailure.build.succeeded);
            testCase.verifyEqual(installFailure.build.reason.code,"install-failed");
            testCase.verifyTrue(contains(installFailure.build.reason.identifier,"InstallFailure"));
            testCase.verifyEqual(arrayfun(@TestFFTWBackend.fileHash,modulePaths),baselineHashes);
            testCase.verifyTrue(FFTWBackend.capabilities().isComplete);

            context.failureInjection = "final-validation-failed";
            validationFailure = TestableFFTWBackend.buildUsing(context);
            testCase.verifyFalse(validationFailure.build.succeeded);
            testCase.verifyEqual(validationFailure.build.reason.code,"final-validation-failed");
            testCase.verifyTrue(contains(validationFailure.build.reason.identifier,"FinalValidationFailure"));
            testCase.verifyEqual(arrayfun(@TestFFTWBackend.fileHash,modulePaths),baselineHashes);
            testCase.verifyTrue(FFTWBackend.capabilities().isComplete);
        end

        function testBuildFailureIsReturnedWithoutWarning(testCase)
            context = TestableFFTWBackend.context();
            context.failureInjection = "compiler-unavailable";
            testCase.verifyWarningFree(@() TestableFFTWBackend.buildUsing(context));
            result = TestableFFTWBackend.buildUsing(context);
            testCase.verifyFalse(result.build.succeeded);
            testCase.verifyEqual(result.build.reason.code,"compiler-unavailable");
            testCase.verifyTrue(contains(result.build.reason.identifier,"CompilerUnavailable"));
        end
    end

    methods (Static, Access=private)
        function value = isEligible(record,batchCount)
            value = record.eligible && any(arrayfun(@(interval) batchCount >= interval.minimumBatchCount && batchCount <= interval.maximumBatchCount,record.intervals));
        end

        function path = repositoryRoot()
            testDirectory = fileparts(mfilename('fullpath'));
            fftwDirectory = fileparts(testDirectory);
            path = fileparts(fileparts(fileparts(fftwDirectory)));
        end

        function paths = modulePaths(context)
            paths = [fullfile(context.sourceDirectory,"fftw_r2c."+context.mexExtension) fullfile(context.sourceDirectory,"fftw_r2r."+context.mexExtension)];
        end

        function hash = fileHash(path)
            fileId = fopen(path,'r');
            cleanup = onCleanup(@() fclose(fileId));
            bytes = fread(fileId,Inf,'*uint8');
            digest = java.security.MessageDigest.getInstance('SHA-256');
            digest.update(bytes);
            hashBytes = typecast(digest.digest(),'uint8');
            hash = string(lower(reshape(dec2hex(hashBytes,2).',1,[])));
            clear cleanup
        end

        function wisdom = canonicalWisdom(wisdom)
            lines = strip(splitlines(string(wisdom)));
            wisdom = sort(lines(strlength(lines) > 0));
        end
    end
end
