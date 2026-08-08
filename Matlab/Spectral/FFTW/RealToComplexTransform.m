classdef RealToComplexTransform < handle
    % Apply lean FFTW real-to-complex and complex-to-real transforms.
    %
    % The final entry of `dims` is stored as a canonical interleaved half
    % spectrum. For example, `dims=[2 1]` produces half-x storage and
    % `dims=[1 2]` produces half-y storage. Inverse transforms use FFTW's
    % unnormalized convention; multiply their result by `scaleFactor` to
    % reproduce MATLAB's normalized inverse.
    %
    % Preserving inverse methods copy the half spectrum exactly once because
    % FFTW destroys multidimensional c2r inputs. The destructive method avoids
    % that copy for uniquely owned arrays and returns the destroyed spectrum,
    % which callers must reassign.
    %
    % ```matlab
    % transform = RealToComplexTransform([128 128 32],dims=[2 1]);
    % spectrum = transform.transformForward(x);
    % xRoundTrip = transform.scaleFactor*transform.transformBack(spectrum);
    % ```
    %
    % - Topic: Create a transform
    % - Topic: Inspect transform properties
    % - Topic: Apply transforms
    % - Topic: Build the backend
    % - Declaration: classdef RealToComplexTransform < handle

    properties (SetAccess=private)
        % Dimensions of the real-valued input array.
        %
        % - Topic: Inspect transform properties
        realSize (1,:) double

        % Dimensions of the canonical interleaved half spectrum.
        %
        % - Topic: Inspect transform properties
        complexSize (1,:) double

        % Factor that normalizes an FFTW complex-to-real result.
        %
        % `scaleFactor` equals the reciprocal product of all transformed
        % dimension lengths.
        %
        % - Topic: Inspect transform properties
        scaleFactor (1,1) double = 1

        % Ordered dimensions transformed by FFTW.
        %
        % The final dimension is the compressed half-spectrum dimension.
        %
        % - Topic: Inspect transform properties
        transformDimensions (1,:) double

        % FFTW new-array alignment policy used by the plan.
        %
        % - Topic: Inspect transform properties
        alignmentMode (1,1) string
    end

    properties (Access=private)
        plan (1,1) uint64 = uint64(0)
    end

    methods
        function self = RealToComplexTransform(sz,options)
            % Create an FFTW transform for a fixed real-array shape.
            %
            % - Topic: Create a transform
            % - Parameter sz: Positive integer dimensions of the real array.
            % - Parameter dims: Nonempty ordered list of distinct transform dimensions.
            % - Parameter planner: FFTW planning mode.
            % - Parameter nCores: Positive FFTW thread count.
            % - Parameter alignmentMode: `"unaligned"` for arbitrary arrays or `"matched"` for strict new-array alignment.
            % - Parameter plannerTimeLimitSeconds: Positive planning limit applied separately to forward and inverse plans.
            % - Returns self: Configured transform object.
            arguments
                sz (1,:) double {mustBeInteger,mustBePositive}
                options.dims (1,:) double {mustBeInteger,mustBePositive} = 1
                options.planner (1,1) string {mustBeMember(options.planner,["estimate","measure","patient","exhaustive"])} = "measure"
                options.nCores (1,1) double {mustBeInteger,mustBePositive} = maxNumCompThreads
                options.alignmentMode (1,1) string {mustBeMember(options.alignmentMode,["matched","unaligned"])} = "unaligned"
                options.plannerTimeLimitSeconds (1,1) double {mustBePositive,mustBeFinite} = 10
            end
            if isempty(options.dims) || numel(unique(options.dims)) ~= numel(options.dims) || any(options.dims > numel(sz))
                error('RealToComplexTransform:InvalidTransformDimensions','dims must be a nonempty ordered list of distinct dimensions within sz.');
            end
            if any(sz(options.dims) == 1)
                error('RealToComplexTransform:SingletonTransformDimension','Transform dimensions must have length greater than one.');
            end
            if exist('fftw_r2c','file') ~= 3
                error('RealToComplexTransform:MexUnavailable','The fftw_r2c MEX backend is unavailable. Run RealToComplexTransform.makeMexFiles first.');
            end

            self.realSize = sz;
            self.transformDimensions = options.dims;
            self.alignmentMode = options.alignmentMode;
            plannerFlags = RealToComplexTransform.plannerFlags(options.planner);
            [self.plan,self.complexSize,self.scaleFactor] = fftw_r2c('create',sz,options.dims,options.nCores,plannerFlags,char(options.alignmentMode),options.plannerTimeLimitSeconds);
        end

        function xbar = transformForward(self,x)
            % Transform a real array into a newly allocated half spectrum.
            %
            % The MEX backend allocates MATLAB-managed storage and transfers it
            % to the returned array without copying the spectrum.
            %
            % - Topic: Apply transforms
            % - Parameter x: Real double array with shape `realSize`.
            % - Returns xbar: Complex half spectrum with shape `complexSize`.
            xbar = fftw_r2c('forward',self.plan,x,'allocating');
        end

        function fbar = transformForwardIntoArray(self,f,fbar)
            % Transform into a caller-provided complex half spectrum.
            %
            % Reassign the returned `fbar`. An aliased output may detach through
            % MATLAB copy-on-write before FFTW executes.
            %
            % - Topic: Apply transforms
            % - Parameter f: Real double array with shape `realSize`.
            % - Parameter fbar: Complex double array with shape `complexSize`.
            % - Returns fbar: Reassigned transformed half spectrum.
            fbar = fftw_r2c('forward',self.plan,f,'preallocated',fbar);
        end

        function x = transformBack(self,xbar)
            % Apply an input-preserving inverse into a new real array.
            %
            % This method performs exactly one half-spectrum copy. Multiply the
            % returned array by `scaleFactor` for MATLAB-style normalization.
            %
            % - Topic: Apply transforms
            % - Parameter xbar: Complex half spectrum with shape `complexSize`.
            % - Returns x: Unnormalized real array with shape `realSize`.
            x = fftw_r2c('inversePreserving',self.plan,xbar,'allocating');
        end

        function x = transformBackIntoArray(self,xbar,x)
            % Apply an input-preserving inverse into a caller-provided array.
            %
            % Reassign the returned `x`. The half spectrum is copied exactly
            % once and remains unchanged.
            %
            % - Topic: Apply transforms
            % - Parameter xbar: Complex half spectrum with shape `complexSize`.
            % - Parameter x: Real double output array with shape `realSize`.
            % - Returns x: Reassigned unnormalized real output.
            x = fftw_r2c('inversePreserving',self.plan,xbar,'preallocated',x);
        end

        function [xbar,x] = transformBackIntoArrayDestructive(self,xbar,x)
            % Apply a destructive inverse without an explicit spectrum copy.
            %
            % Both outputs must be reassigned. Uniquely owned inputs retain
            % their pointers; aliased inputs may detach through MATLAB
            % copy-on-write. Multiply `x` by `scaleFactor` to normalize it.
            %
            % - Topic: Apply transforms
            % - Parameter xbar: Complex half spectrum that FFTW may destroy.
            % - Parameter x: Real double output array with shape `realSize`.
            % - Returns xbar: Reassigned destroyed spectrum.
            % - Returns x: Reassigned unnormalized real output.
            [xbar,x] = fftw_r2c('inverseDestructive',self.plan,xbar,x);
        end

        function delete(self)
            if self.plan ~= 0 && exist('fftw_r2c','file') == 3
                fftw_r2c('free',self.plan);
                self.plan = uint64(0);
            end
        end
    end

    methods (Static)
        function makeMexFiles(fftwlibpath)
            % Build the hardened backend against an FFTW library.
            %
            % The default is the active MATLAB installation's bundled FFTW.
            % Source and output paths are resolved relative to this class file,
            % so the build is independent of the current working directory.
            %
            % - Topic: Build the backend
            % - Parameter fftwlibpath: Path to an FFTW-compatible dynamic library.
            arguments
                fftwlibpath (1,1) string = fullfile(matlabroot,'bin',computer('arch'),'libmwfftw3.3.dylib')
            end
            sourceDirectory = fileparts(mfilename('fullpath'));
            sourcePath = fullfile(sourceDirectory,'fftw_r2c.cpp');
            includeArgument = "-I" + sourceDirectory;
            clear fftw_r2c
            mex('-R2018a','-outdir',sourceDirectory,'-output','fftw_r2c',includeArgument,sourcePath,fftwlibpath);
            rehash
        end
    end

    methods (Static, Access=private)
        function flags = plannerFlags(planner)
            switch planner
                case "estimate"
                    flags = bitshift(1,6);
                case "measure"
                    flags = 0;
                case "patient"
                    flags = bitshift(1,5);
                case "exhaustive"
                    flags = bitshift(1,3);
            end
        end
    end
end
