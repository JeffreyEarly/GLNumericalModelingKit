classdef TestRealToRealTransform < matlab.unittest.TestCase
    methods (TestClassSetup)
        function buildBackend(testCase)
            fftwDirectory = TestRealToRealTransform.fftwDirectory;
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fftwDirectory));
            addpath(fullfile(TestRealToRealTransform.repositoryRoot,'Matlab','Spectral'));
            RealToRealTransform.makeMexFiles;
        end
    end

    methods (Test)
        function testMatrixAgreementAndRoundTrips(testCase)
            configurations = {
                [9 4], 1, "cosine", "real"
                [9 4], 1, "sine", "real"
                [9 4], 1, "cosine", "complex"
                [9 4], 1, "sine", "complex"
                [8 3], 1, "cosine", "real"
                [8 3], 1, "sine", "real"
                [3 9 4], 2, "cosine", "complex"
                [3 9 4], 2, "sine", "complex"
                [9 4 1], 1, "cosine", "real"
                };

            for iConfiguration = 1:size(configurations,1)
                sz = configurations{iConfiguration,1};
                dimension = configurations{iConfiguration,2};
                transformType = configurations{iConfiguration,3};
                dataType = configurations{iConfiguration,4};
                rng(iConfiguration,'twister');
                values = randn(sz);
                if dataType == "complex"
                    values = complex(values,randn(sz));
                end
                if transformType == "sine"
                    if dataType == "complex"
                        values = TestRealToRealTransform.setEndpoints(values,dimension,3-2i,5+4i);
                    else
                        values = TestRealToRealTransform.setEndpoints(values,dimension,3,5);
                    end
                end

                transform = RealToRealTransform(sz,dims=dimension,transform=transformType,dataType=dataType,planner="estimate",nCores=1);
                cleanup = onCleanup(@() delete(transform));
                [forwardMatrix,backMatrix] = TestRealToRealTransform.referenceMatrices(sz(dimension),transformType);
                expected = TestRealToRealTransform.applyMatrix(values,forwardMatrix,dimension);
                coefficients = transform.transformForward(values);
                testCase.verifyEqual(size(coefficients),TestRealToRealTransform.canonicalSize(transform.spectralSize));
                testCase.verifyLessThanOrEqual(TestRealToRealTransform.relativeError(coefficients,expected),1e-12);

                preallocated = TestRealToRealTransform.zerosLike(transform.spectralSize,dataType);
                preallocated = transform.transformForwardIntoArray(values,preallocated);
                testCase.verifyLessThanOrEqual(TestRealToRealTransform.relativeError(preallocated,coefficients),1e-12);

                expectedBack = TestRealToRealTransform.applyMatrix(coefficients,backMatrix,dimension);
                roundTrip = transform.transformBack(coefficients);
                testCase.verifyLessThanOrEqual(TestRealToRealTransform.relativeError(roundTrip,expectedBack),1e-12);
                output = TestRealToRealTransform.zerosLike(sz,dataType);
                output = transform.transformBackIntoArray(coefficients,output);
                testCase.verifyLessThanOrEqual(TestRealToRealTransform.relativeError(output,roundTrip),1e-12);
                if transformType == "sine"
                    zeroEndpointValues = TestRealToRealTransform.setEndpoints(values,dimension,0,0);
                    testCase.verifyLessThanOrEqual(TestRealToRealTransform.relativeError(roundTrip,zeroEndpointValues),1e-12);
                    endpoints = TestRealToRealTransform.endpointValues(roundTrip,dimension);
                    testCase.verifyEqual(endpoints,zeros(size(endpoints),'like',endpoints));
                else
                    testCase.verifyLessThanOrEqual(TestRealToRealTransform.relativeError(roundTrip,values),1e-12);
                end
                testCase.verifyEqual(transform.scaleFactor,1);
                clear cleanup
            end
        end

        function testOwnershipInstrumentationAndLifetime(testCase)
            fftw_r2r('resetLifetime');
            sz = [9 4];
            [plan,spectralSize] = fftw_r2r('create',sz,1,'cosine','complex',1,64,'unaligned',1);
            cleanup = onCleanup(@() fftw_r2r('free',plan));
            info = fftw_r2r('planInfo',plan);
            testCase.verifyEqual(info(1:3),[prod(sz) prod(spectralSize) 0]);
            input = complex(randn(sz),randn(sz));

            output = fftw_r2r('forward',plan,input,'allocating');
            [metrics,pointers] = fftw_r2r('metrics',plan);
            testCase.verifyEqual(metrics(9:10),[0 0]);
            testCase.verifyEqual(pointers(3),pointers(4));
            testCase.verifyEqual(fftw_r2r('pointer',output),pointers(4));

            [uniqueMetrics,uniquePointers,preallocated] = TestRealToRealTransform.uniquePreallocatedCall(plan,input,spectralSize);
            testCase.verifyEqual(uniqueMetrics(9:10),[0 0]);
            testCase.verifyEqual(uniquePointers(2),uniquePointers(3));
            testCase.verifyEqual(fftw_r2r('pointer',preallocated),uniquePointers(3));
            testCase.verifyLessThanOrEqual(TestRealToRealTransform.relativeError(preallocated,output),1e-12);

            [aliasMetrics,alias,detached] = TestRealToRealTransform.aliasedPreallocatedCall(plan,input,spectralSize);
            testCase.verifyEqual(aliasMetrics(9),1);
            testCase.verifyEqual(aliasMetrics(10),16*prod(spectralSize));
            testCase.verifyEqual(alias,complex(zeros(spectralSize)));
            testCase.verifyNotEqual(detached,alias);

            fftw_r2r('free',plan);
            clear cleanup
            fftw_r2r('free',plan);
            lifetime = fftw_r2r('lifetime');
            testCase.verifyEqual(lifetime(1:3),[1 1 0]);
            testCase.verifyEqual(lifetime(4),lifetime(5));
            testCase.verifyEqual(lifetime(6),0);
        end

        function testCompatibilityAndAlignment(testCase)
            input = randn(9,4);
            production = RealToRealTransform([9 4],dims=1,transform="cosine",planner="estimate",nCores=1);
            productionCleanup = onCleanup(@() delete(production));
            compatibility = RealToRealTransformMexFFTW([9 4],dim=1,transform="cosine",planner="estimate",nCores=1);
            compatibilityCleanup = onCleanup(@() delete(compatibility));
            testCase.verifyLessThanOrEqual(TestRealToRealTransform.relativeError(production.transformForward(input),compatibility.transformForward(input)),1e-12);
            clear compatibilityCleanup productionCleanup

            [matchedAccepted,mismatchRejected,unalignedAccepted] = fftw_r2r('alignmentSelfTest');
            testCase.verifyTrue(matchedAccepted);
            testCase.verifyTrue(mismatchRejected);
            testCase.verifyTrue(unalignedAccepted);
            matched = RealToRealTransform([9 4],dims=1,transform="cosine",planner="estimate",nCores=1,alignmentMode="matched");
            matchedCleanup = onCleanup(@() delete(matched));
            testCase.verifyLessThanOrEqual(TestRealToRealTransform.relativeError(matched.transformBack(matched.transformForward(input)),input),1e-12);
            clear matchedCleanup
        end

        function testErrorsCleanupAndBuildPaths(testCase)
            testCase.verifyError(@() RealToRealTransform([2 4],dims=1,transform="cosine"),'RealToRealTransform:TransformTooShort');
            testCase.verifyError(@() RealToRealTransform([9 4],dims=3,transform="cosine"),'RealToRealTransform:InvalidTransformDimensions');

            fftw_r2r('resetLifetime');
            transform = RealToRealTransform([9 4],dims=1,transform="sine",planner="estimate",nCores=1);
            testCase.verifyError(@() transform.transformForward(zeros(8,4)),'RealToRealTransform:DimensionMismatch');
            testCase.verifyError(@() transform.transformForward(complex(zeros(9,4))),'RealToRealTransform:DataTypeMismatch');
            valid = transform.transformForward(zeros(9,4));
            testCase.verifySize(valid,[7 4]);
            delete(transform);
            lifetime = fftw_r2r('lifetime');
            testCase.verifyEqual(lifetime(1:3),[1 1 0]);

            import matlab.unittest.fixtures.TemporaryFolderFixture
            fixture = testCase.applyFixture(TemporaryFolderFixture);
            previousDirectory = pwd;
            directoryCleanup = onCleanup(@() cd(previousDirectory));
            cd(TestRealToRealTransform.repositoryRoot);
            RealToRealTransform.makeMexFiles;
            cd(fixture.Folder);
            RealToRealTransformMexFFTW.makeMexFiles;
            testCase.verifyEqual(exist('fftw_r2r','file'),3);
            clear directoryCleanup
        end
    end

    methods (Static, Access=private)
        function [forwardMatrix,backMatrix] = referenceMatrices(n,transformType)
            if transformType == "cosine"
                forwardMatrix = CosineTransformForwardMatrix(n);
                backMatrix = CosineTransformBackMatrix(n);
            else
                forwardMatrix = SineTransformForwardMatrix(n);
                backMatrix = SineTransformBackMatrix(n);
            end
        end

        function output = applyMatrix(input,matrix,dimension)
            order = [dimension 1:(dimension-1) (dimension+1):ndims(input)];
            permuted = permute(input,order);
            outputSize = size(permuted);
            transformed = matrix*reshape(permuted,outputSize(1),[]);
            outputSize(1) = size(matrix,1);
            output = ipermute(reshape(transformed,outputSize),order);
        end

        function values = setEndpoints(values,dimension,leftValue,rightValue)
            left = repmat({':'},1,ndims(values));
            right = left;
            left{dimension} = 1;
            right{dimension} = size(values,dimension);
            values(left{:}) = leftValue;
            values(right{:}) = rightValue;
        end

        function values = endpointValues(input,dimension)
            left = repmat({':'},1,ndims(input));
            right = left;
            left{dimension} = 1;
            right{dimension} = size(input,dimension);
            values = [input(left{:}); input(right{:})];
        end

        function output = zerosLike(sz,dataType)
            output = zeros(sz);
            if dataType == "complex"
                output = complex(output);
            end
        end

        function sz = canonicalSize(sz)
            while numel(sz) > 2 && sz(end) == 1
                sz(end) = [];
            end
        end

        function error = relativeError(actual,expected)
            error = norm(actual(:)-expected(:),inf)/max(norm(expected(:),inf),eps);
        end

        function [metrics,pointers,output] = uniquePreallocatedCall(plan,input,outputSize)
            output = complex(zeros(outputSize));
            output = fftw_r2r('forward',plan,input,'preallocated',output);
            [metrics,pointers] = fftw_r2r('metrics',plan);
        end

        function [metrics,alias,output] = aliasedPreallocatedCall(plan,input,outputSize)
            output = complex(zeros(outputSize));
            alias = output;
            output = fftw_r2r('forward',plan,input,'preallocated',output);
            metrics = fftw_r2r('metrics',plan);
        end

        function path = fftwDirectory()
            path = fileparts(fileparts(mfilename('fullpath')));
        end

        function path = repositoryRoot()
            path = fileparts(fileparts(fileparts(TestRealToRealTransform.fftwDirectory)));
        end
    end
end
