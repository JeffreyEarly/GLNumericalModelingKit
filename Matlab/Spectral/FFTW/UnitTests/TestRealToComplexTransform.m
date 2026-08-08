classdef TestRealToComplexTransform < matlab.unittest.TestCase
    methods (TestClassSetup)
        function buildBackend(testCase)
            fftwDirectory = TestRealToComplexTransform.fftwDirectory;
            addpathCleanup = testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fftwDirectory)); %#ok<NASGU>
            RealToComplexTransform.makeMexFiles;
        end
    end

    methods (Test)
        function testOrderedDimensionsLayoutsAndNormalization(testCase)
            cases = {
                [17 5], 1, [9 5]
                [16 8 4], [1 2], [16 5 4]
                [16 8 4], [2 1], [9 8 4]
                [16 8 4], [3 1], [9 8 4]
                [8 6 4], [2 3 1], [5 6 4]
                [16 8 1], 1, [9 8 1]
                };

            for iCase = 1:size(cases,1)
                sz = cases{iCase,1};
                dimensions = cases{iCase,2};
                expectedComplexSize = cases{iCase,3};
                transform = RealToComplexTransform(sz,dims=dimensions,planner="estimate",nCores=1);
                cleanup = onCleanup(@() delete(transform));
                testCase.verifyEqual(transform.complexSize,expectedComplexSize);
                testCase.verifyEqual(transform.scaleFactor,1/prod(sz(dimensions)),'AbsTol',eps);

                rng(iCase,'twister');
                input = randn(sz);
                expectedSpectrum = TestRealToComplexTransform.matlabHalfSpectrum(input,dimensions);
                spectrum = transform.transformForward(input);
                testCase.verifyLessThanOrEqual(TestRealToComplexTransform.relativeError(spectrum,expectedSpectrum),1e-12);

                preallocatedSpectrum = complex(zeros(expectedComplexSize));
                preallocatedSpectrum = transform.transformForwardIntoArray(input,preallocatedSpectrum);
                testCase.verifyLessThanOrEqual(TestRealToComplexTransform.relativeError(preallocatedSpectrum,spectrum),1e-12);

                preservedSpectrum = complex(real(spectrum),imag(spectrum));
                spectrumBefore = preservedSpectrum;
                output = transform.transformBack(preservedSpectrum);
                testCase.verifyEqual(preservedSpectrum,spectrumBefore);
                testCase.verifyLessThanOrEqual(TestRealToComplexTransform.relativeError(transform.scaleFactor*output,input),1e-12);

                preallocatedOutput = zeros(sz);
                preallocatedOutput = transform.transformBackIntoArray(preservedSpectrum,preallocatedOutput);
                testCase.verifyEqual(preservedSpectrum,spectrumBefore);
                testCase.verifyLessThanOrEqual(TestRealToComplexTransform.relativeError(transform.scaleFactor*preallocatedOutput,input),1e-12);
                clear cleanup
            end
            testCase.verifyFalse(any(strcmp(properties('RealToComplexTransform'),'scratch')));
        end

        function testOwnershipCopyAndLifetimeContracts(testCase)
            fftw_r2c('resetLifetime');
            sz = [16 8 4];
            dimensions = [2 1];
            [plan,complexSize] = fftw_r2c('create',sz,dimensions,1,64,'unaligned',1);
            cleanup = onCleanup(@() fftw_r2c('free',plan));
            info = fftw_r2c('planInfo',plan);
            testCase.verifyEqual(info(1),prod(sz));
            testCase.verifyEqual(info(2),prod(complexSize));
            testCase.verifyEqual(info(3),16*prod(complexSize));

            rng(41,'twister');
            input = randn(sz);
            spectrum = fftw_r2c('forward',plan,input,'allocating');
            [metrics,pointers] = fftw_r2c('metrics',plan);
            returnedPointer = fftw_r2c('pointer',spectrum);
            testCase.verifyEqual(metrics(9:12),zeros(1,4));
            testCase.verifyEqual(pointers(4),pointers(5));
            testCase.verifyEqual(returnedPointer,pointers(5));

            [preallocatedMetrics,preallocatedPointers,preallocatedSpectrum] = TestRealToComplexTransform.preallocatedForwardCall(plan,input,complexSize);
            testCase.verifyEqual(preallocatedMetrics(9:12),zeros(1,4));
            testCase.verifyEqual(preallocatedPointers(3),preallocatedPointers(4));
            testCase.verifyEqual(fftw_r2c('pointer',preallocatedSpectrum),preallocatedPointers(4));
            testCase.verifyLessThanOrEqual(TestRealToComplexTransform.relativeError(preallocatedSpectrum,spectrum),1e-12);

            preservingInput = complex(real(spectrum),imag(spectrum));
            preservingBefore = preservingInput;
            output = fftw_r2c('inversePreserving',plan,preservingInput,'allocating'); %#ok<NASGU>
            metrics = fftw_r2c('metrics',plan);
            testCase.verifyEqual(metrics(9),1);
            testCase.verifyEqual(metrics(10),16*prod(complexSize));
            testCase.verifyEqual(preservingInput,preservingBefore);

            [preallocatedInverseMetrics,preallocatedInversePointers,preservedOutput] = TestRealToComplexTransform.preallocatedPreservingCall(plan,preservingInput,sz);
            testCase.verifyEqual(preallocatedInverseMetrics(9),1);
            testCase.verifyEqual(preallocatedInverseMetrics(10),16*prod(complexSize));
            testCase.verifyEqual(preallocatedInverseMetrics(11:12),[0 0]);
            testCase.verifyEqual(preallocatedInversePointers(3),preallocatedInversePointers(4));
            testCase.verifyLessThanOrEqual(TestRealToComplexTransform.relativeError(preservedOutput/prod(sz(dimensions)),input),1e-12);

            [destructiveMetrics,destructivePointers,destroyed,realOutput] = TestRealToComplexTransform.uniqueDestructiveCall(plan,spectrum,sz);
            testCase.verifyEqual(destructiveMetrics(9:12),zeros(1,4));
            testCase.verifyEqual(destructivePointers(1),destructivePointers(2));
            testCase.verifyEqual(fftw_r2c('pointer',destroyed),destructivePointers(2));
            testCase.verifyLessThanOrEqual(TestRealToComplexTransform.relativeError(realOutput/prod(sz(dimensions)),input),1e-12);

            [aliasedMetrics,spectrumAlias,outputAlias,destroyedAliased,outputAliased] = TestRealToComplexTransform.aliasedDestructiveCall(plan,spectrum,sz);
            testCase.verifyEqual(aliasedMetrics(11),2);
            testCase.verifyEqual(aliasedMetrics(12),16*prod(complexSize)+8*prod(sz));
            testCase.verifyEqual(spectrumAlias,spectrum);
            testCase.verifyEqual(outputAlias,zeros(sz));
            testCase.verifyNotEqual(destroyedAliased,spectrumAlias);
            testCase.verifyNotEqual(outputAliased,outputAlias);

            fftw_r2c('free',plan);
            clear cleanup
            fftw_r2c('free',plan);
            lifetime = fftw_r2c('lifetime');
            testCase.verifyEqual(lifetime(1:3),[1 1 0]);
            testCase.verifyEqual(lifetime(4),lifetime(5));
            testCase.verifyEqual(lifetime(6),0);
        end

        function testAlignmentModes(testCase)
            [matchedAccepted,mismatchRejected,unalignedAccepted] = fftw_r2c('alignmentSelfTest');
            testCase.verifyTrue(matchedAccepted);
            testCase.verifyTrue(mismatchRejected);
            testCase.verifyTrue(unalignedAccepted);

            input = randn(16,8,4);
            unaligned = RealToComplexTransform([16 8 4],dims=[2 1],planner="estimate",nCores=1);
            unalignedCleanup = onCleanup(@() delete(unaligned));
            unalignedSpectrum = unaligned.transformForward(input);
            testCase.verifyLessThanOrEqual(TestRealToComplexTransform.relativeError(unalignedSpectrum,TestRealToComplexTransform.matlabHalfSpectrum(input,[2 1])),1e-12);
            clear unalignedCleanup

            matched = RealToComplexTransform([16 8 4],dims=[2 1],planner="estimate",nCores=1,alignmentMode="matched");
            matchedCleanup = onCleanup(@() delete(matched));
            matchedSpectrum = matched.transformForward(input);
            testCase.verifyLessThanOrEqual(TestRealToComplexTransform.relativeError(matchedSpectrum,TestRealToComplexTransform.matlabHalfSpectrum(input,[2 1])),1e-12);
            clear matchedCleanup
        end

        function testClassDestructiveReassignmentAndAliases(testCase)
            sz = [16 8 4];
            transform = RealToComplexTransform(sz,dims=[2 1],planner="estimate",nCores=1);
            cleanup = onCleanup(@() delete(transform));
            input = randn(sz);
            spectrum = transform.transformForward(input);

            [beforeSpectrum,beforeOutput,afterSpectrum,afterOutput,destroyed,output] = TestRealToComplexTransform.uniqueClassDestructiveCall(transform,spectrum,sz);
            testCase.verifyEqual(afterSpectrum,beforeSpectrum);
            testCase.verifyEqual(afterOutput,beforeOutput);
            testCase.verifyEqual(fftw_r2c('pointer',destroyed),beforeSpectrum);
            testCase.verifyEqual(fftw_r2c('pointer',output),beforeOutput);
            testCase.verifyLessThanOrEqual(TestRealToComplexTransform.relativeError(transform.scaleFactor*output,input),1e-12);

            [spectrumAlias,outputAlias,destroyedAliased,outputAliased] = TestRealToComplexTransform.aliasedClassDestructiveCall(transform,spectrum,sz);
            testCase.verifyEqual(spectrumAlias,spectrum);
            testCase.verifyEqual(outputAlias,zeros(sz));
            testCase.verifyNotEqual(destroyedAliased,spectrumAlias);
            testCase.verifyNotEqual(outputAliased,outputAlias);
            clear cleanup
        end

        function testErrorsAndCleanup(testCase)
            testCase.verifyError(@() RealToComplexTransform([8 8],dims=[1 1]),'RealToComplexTransform:InvalidTransformDimensions');
            testCase.verifyError(@() RealToComplexTransform([8 8],dims=3),'RealToComplexTransform:InvalidTransformDimensions');
            testCase.verifyError(@() RealToComplexTransform([8 1],dims=2),'RealToComplexTransform:SingletonTransformDimension');

            fftw_r2c('resetLifetime');
            for iTransform = 1:3
                transform = RealToComplexTransform([8 8],dims=[2 1],planner="estimate",nCores=1);
                testCase.verifyError(@() transform.transformForward(zeros(7,8)),'RealToComplexTransform:DimensionMismatch');
                testCase.verifyError(@() transform.transformForward(complex(zeros(8,8))),'RealToComplexTransform:InvalidRealInput');
                valid = transform.transformForward(zeros(8,8));
                testCase.verifySize(valid,[5 8]);
                delete(transform);
            end
            lifetime = fftw_r2c('lifetime');
            testCase.verifyEqual(lifetime(1:3),[3 3 0]);
            testCase.verifyEqual(lifetime(6),0);
        end

        function testBuildIsIndependentOfCurrentDirectory(testCase)
            import matlab.unittest.fixtures.TemporaryFolderFixture

            fixture = testCase.applyFixture(TemporaryFolderFixture);
            previousDirectory = pwd;
            directoryCleanup = onCleanup(@() cd(previousDirectory));
            cd(TestRealToComplexTransform.repositoryRoot);
            RealToComplexTransform.makeMexFiles;
            cd(fixture.Folder);
            RealToComplexTransform.makeMexFiles;
            testCase.verifyEqual(exist('fftw_r2c','file'),3);
            clear directoryCleanup
        end
    end

    methods (Static, Access=private)
        function spectrum = matlabHalfSpectrum(input,dimensions)
            spectrum = input;
            for dimension = dimensions
                spectrum = fft(spectrum,[],dimension);
            end
            indices = repmat({':'},1,ndims(spectrum));
            compressedDimension = dimensions(end);
            indices{compressedDimension} = 1:(floor(size(input,compressedDimension)/2)+1);
            spectrum = spectrum(indices{:});
        end

        function error = relativeError(actual,expected)
            error = norm(actual(:)-expected(:),inf)/max(norm(expected(:),inf),eps);
        end

        function [metrics,pointers,spectrum,output] = uniqueDestructiveCall(plan,referenceSpectrum,sz)
            spectrum = complex(real(referenceSpectrum),imag(referenceSpectrum));
            output = zeros(sz);
            [spectrum,output] = fftw_r2c('inverseDestructive',plan,spectrum,output);
            [metrics,pointers] = fftw_r2c('metrics',plan);
        end

        function [metrics,pointers,spectrum] = preallocatedForwardCall(plan,input,complexSize)
            spectrum = complex(zeros(complexSize));
            spectrum = fftw_r2c('forward',plan,input,'preallocated',spectrum);
            [metrics,pointers] = fftw_r2c('metrics',plan);
        end

        function [metrics,pointers,output] = preallocatedPreservingCall(plan,spectrum,sz)
            output = zeros(sz);
            output = fftw_r2c('inversePreserving',plan,spectrum,'preallocated',output);
            [metrics,pointers] = fftw_r2c('metrics',plan);
        end

        function [metrics,spectrumAlias,outputAlias,spectrum,output] = aliasedDestructiveCall(plan,referenceSpectrum,sz)
            spectrum = complex(real(referenceSpectrum),imag(referenceSpectrum));
            output = zeros(sz);
            spectrumAlias = spectrum;
            outputAlias = output;
            [spectrum,output] = fftw_r2c('inverseDestructive',plan,spectrum,output);
            metrics = fftw_r2c('metrics',plan);
        end

        function [beforeSpectrum,beforeOutput,afterSpectrum,afterOutput,spectrum,output] = uniqueClassDestructiveCall(transform,referenceSpectrum,sz)
            spectrum = complex(real(referenceSpectrum),imag(referenceSpectrum));
            output = zeros(sz);
            beforeSpectrum = fftw_r2c('pointer',spectrum);
            beforeOutput = fftw_r2c('pointer',output);
            [spectrum,output] = transform.transformBackIntoArrayDestructive(spectrum,output);
            afterSpectrum = fftw_r2c('pointer',spectrum);
            afterOutput = fftw_r2c('pointer',output);
        end

        function [spectrumAlias,outputAlias,spectrum,output] = aliasedClassDestructiveCall(transform,referenceSpectrum,sz)
            spectrum = complex(real(referenceSpectrum),imag(referenceSpectrum));
            output = zeros(sz);
            spectrumAlias = spectrum;
            outputAlias = output;
            [spectrum,output] = transform.transformBackIntoArrayDestructive(spectrum,output);
        end

        function path = fftwDirectory()
            path = fileparts(fileparts(mfilename('fullpath')));
        end

        function path = repositoryRoot()
            path = fileparts(fileparts(fileparts(TestRealToComplexTransform.fftwDirectory)));
        end
    end
end
