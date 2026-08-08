classdef RealToRealTransform < handle
    % Apply normalized FFTW DCT-I and DST-I transforms.
    %
    % RealToRealTransform follows the same ownership and build conventions
    % as RealToComplexTransform. The configured dims entry identifies the
    % single transformed dimension. Cosine transforms retain its length;
    % sine transforms omit the two physical endpoint values spectrally.
    %
    % Both real and interleaved-complex arrays are supported. Complex plans
    % transform real and imaginary components as FFTW batches without
    % constructing separate MATLAB arrays.
    %
    % Example:
    %   transform = RealToRealTransform([65 128],dims=1,transform="cosine",dataType="complex");
    %   coefficients = transform.transformForward(values);
    %   roundTrip = transform.transformBack(coefficients);
    %
    % - Topic: Create a transform
    % - Topic: Inspect transform properties
    % - Topic: Apply transforms
    % - Topic: Build the backend
    % - Declaration: classdef RealToRealTransform < handle

    properties (SetAccess=private)
        % Dimensions of the physical input array.
        %
        % - Topic: Inspect transform properties
        realSize (1,:) double

        % Dimensions of the spectral coefficient array.
        %
        % A sine transform shortens the transformed dimension by two.
        %
        % - Topic: Inspect transform properties
        spectralSize (1,:) double

        % Additional caller normalization required after a transform.
        %
        % Results already follow GL conventions, so this value is one.
        %
        % - Topic: Inspect transform properties
        scaleFactor (1,1) double = 1

        % Single dimension transformed by FFTW.
        %
        % - Topic: Inspect transform properties
        transformDimensions (1,1) double

        % Real-to-real transform family, "cosine" or "sine".
        %
        % - Topic: Inspect transform properties
        transformType (1,1) string

        % Array storage type, "real" or "complex".
        %
        % - Topic: Inspect transform properties
        dataType (1,1) string

        % FFTW new-array alignment policy used by the plan.
        %
        % - Topic: Inspect transform properties
        alignmentMode (1,1) string
    end

    properties (Access=private)
        plan (1,1) uint64 = uint64(0)
    end

    methods
        function self = RealToRealTransform(sz,options)
            % Create a normalized DCT-I or DST-I transform.
            %
            % - Topic: Create a transform
            % - Parameter sz: Positive integer dimensions of the physical array.
            % - Parameter dims: Single nonsingleton transform dimension.
            % - Parameter transform: "cosine" for DCT-I or "sine" for DST-I.
            % - Parameter dataType: "real" or "complex" input and output storage.
            % - Parameter planner: FFTW planning mode.
            % - Parameter nCores: Positive FFTW thread count.
            % - Parameter alignmentMode: "unaligned" or strict "matched" execution.
            % - Parameter plannerTimeLimitSeconds: Positive limit for each plan.
            % - Returns self: Configured transform object.
            arguments
                sz (1,:) double {mustBeInteger,mustBePositive}
                options.dims (1,1) double {mustBeInteger,mustBePositive} = 1
                options.transform (1,1) string {mustBeMember(options.transform,["cosine","sine"])}
                options.dataType (1,1) string {mustBeMember(options.dataType,["real","complex"])} = "real"
                options.planner (1,1) string {mustBeMember(options.planner,["estimate","measure","patient","exhaustive"])} = "measure"
                options.nCores (1,1) double {mustBeInteger,mustBePositive} = maxNumCompThreads
                options.alignmentMode (1,1) string {mustBeMember(options.alignmentMode,["matched","unaligned"])} = "unaligned"
                options.plannerTimeLimitSeconds (1,1) double {mustBePositive,mustBeFinite} = 10
            end
            if options.dims > numel(sz)
                error('RealToRealTransform:InvalidTransformDimensions','dims must identify one dimension within sz.');
            end
            if sz(options.dims) < 3
                error('RealToRealTransform:TransformTooShort','DCT-I/DST-I transforms require at least three physical points.');
            end
            if exist('fftw_r2r','file') ~= 3
                error('RealToRealTransform:MexUnavailable','The fftw_r2r MEX backend is unavailable. Run RealToRealTransform.makeMexFiles first.');
            end

            self.realSize = sz;
            self.spectralSize = sz;
            if options.transform == "sine"
                self.spectralSize(options.dims) = self.spectralSize(options.dims)-2;
            end
            self.transformDimensions = options.dims;
            self.transformType = options.transform;
            self.dataType = options.dataType;
            self.alignmentMode = options.alignmentMode;
            flags = RealToRealTransform.plannerFlags(options.planner);
            [self.plan,mexSpectralSize,self.scaleFactor] = fftw_r2r('create',sz,options.dims,char(options.transform),char(options.dataType),options.nCores,flags,char(options.alignmentMode),options.plannerTimeLimitSeconds);
            if ~isequal(self.spectralSize,mexSpectralSize)
                fftw_r2r('free',self.plan);
                self.plan = uint64(0);
                error('RealToRealTransform:BackendShapeMismatch','The MEX backend returned an unexpected spectral shape.');
            end
        end

        function coefficients = transformForward(self,values)
            % Transform a physical array into newly allocated coefficients.
            %
            % Sine input endpoints are ignored. The returned coefficients are
            % already normalized to match the GL forward-transform matrices.
            %
            % - Topic: Apply transforms
            % - Parameter values: Physical array with shape realSize.
            % - Returns coefficients: Spectral array with shape spectralSize.
            coefficients = fftw_r2r('forward',self.plan,values,'allocating');
        end

        function coefficients = transformForwardIntoArray(self,values,coefficients)
            % Transform into a caller-provided coefficient array.
            %
            % Reassign the returned output. Aliased storage may detach through
            % MATLAB copy-on-write before FFTW executes.
            %
            % - Topic: Apply transforms
            % - Parameter values: Physical array with shape realSize.
            % - Parameter coefficients: Output array with shape spectralSize.
            % - Returns coefficients: Reassigned normalized coefficients.
            coefficients = fftw_r2r('forward',self.plan,values,'preallocated',coefficients);
        end

        function values = transformBack(self,coefficients)
            % Transform coefficients into a newly allocated physical array.
            %
            % A sine inverse restores exact zero physical endpoints.
            %
            % - Topic: Apply transforms
            % - Parameter coefficients: Spectral array with shape spectralSize.
            % - Returns values: Physical array with shape realSize.
            values = fftw_r2r('back',self.plan,coefficients,'allocating');
        end

        function values = transformBackIntoArray(self,coefficients,values)
            % Transform coefficients into a caller-provided physical array.
            %
            % Reassign the returned output. Aliased storage may detach through
            % MATLAB copy-on-write before FFTW executes.
            %
            % - Topic: Apply transforms
            % - Parameter coefficients: Spectral array with shape spectralSize.
            % - Parameter values: Output array with shape realSize.
            % - Returns values: Reassigned normalized physical array.
            values = fftw_r2r('back',self.plan,coefficients,'preallocated',values);
        end

        function delete(self)
            if self.plan ~= 0 && exist('fftw_r2r','file') == 3
                fftw_r2r('free',self.plan);
                self.plan = uint64(0);
            end
        end
    end

    methods (Hidden)
        function metrics = backendMetrics(self)
            % Return private backend instrumentation for tests and benchmarks.
            metrics = fftw_r2r('metrics',self.plan);
        end
    end

    methods (Static)
        function makeMexFiles(fftwlibpath)
            % Build the production backend against an FFTW library.
            %
            % Prefer `FFTWBackend.build()` for release-aware bundled-library
            % discovery, validation, and transactional installation. This
            % method remains available for compatibility and experiments.
            %
            % The default resolves MATLAB's bundled FFTW library. Source and
            % output paths are independent of the current working directory.
            %
            % - Topic: Build the backend
            % - Parameter fftwlibpath: FFTW-compatible dynamic library path.
            arguments
                fftwlibpath (1,1) string = fullfile(matlabroot,'bin',computer('arch'),'libmwfftw3.3.dylib')
            end
            sourceDirectory = fileparts(mfilename('fullpath'));
            sourcePath = fullfile(sourceDirectory,'fftw_r2r.cpp');
            includeArgument = "-I" + sourceDirectory;
            clear fftw_r2r
            mex('-R2018a','-outdir',sourceDirectory,'-output','fftw_r2r',includeArgument,sourcePath,fftwlibpath);
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
