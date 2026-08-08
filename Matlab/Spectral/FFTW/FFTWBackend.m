classdef FFTWBackend
    % Discover, build, and validate the optional bundled-FFTW backend.
    %
    % Capability queries never build the MEX modules and return structured
    % availability rather than throwing for unsupported or failed backends.
    % The build method stages and validates both production modules before
    % installing them beside this class.
    %
    % ```matlab
    % capabilities = FFTWBackend.capabilities();
    % if capabilities.build.isRequired
    %     capabilities = FFTWBackend.build();
    % end
    % ```
    %
    % - Topic: Discover backend support
    % - Topic: Build the backend
    % - Declaration: classdef FFTWBackend

    methods (Static)
        function capabilities = capabilities()
            % Return validated FFTW capabilities without building.
            %
            % Expected availability failures are returned in `reason` and
            % per-module or per-feature reason records.
            %
            % - Topic: Discover backend support
            % - Returns capabilities: Scalar structured capability record.
            context = FFTWBackend.defaultContext;
            capabilities = FFTWBackend.capabilitiesWithContext(context);
        end

        function capabilities = build()
            % Build and validate the bundled-FFTW production modules.
            %
            % The active MATLAB installation is the only library source.
            % Expected build and availability failures are returned rather
            % than thrown.
            %
            % - Topic: Build the backend
            % - Returns capabilities: Scalar structured capability record.
            context = FFTWBackend.defaultContext;
            capabilities = FFTWBackend.buildWithContext(context);
        end
    end

    methods (Static, Access=protected)
        function context = defaultContext()
            % Protected so tests can exercise provider-neutral orchestration.
            sourceDirectory = fileparts(mfilename('fullpath'));
            context.release = string(version('-release'));
            context.architecture = string(computer('arch'));
            context.matlabVersion = string(version);
            context.matlabRoot = string(matlabroot);
            context.sourceDirectory = string(sourceDirectory);
            context.moduleDirectory = string(sourceDirectory);
            context.mexExtension = string(mexext);
            context.provider = FFTWBackend.bundledProvider(context);
            context.errorTolerance = 1e-12;
            context.failureInjection = "";
        end

        function capabilities = capabilitiesWithContext(context)
            state = FFTWBackend.captureMatlabState;
            stateCleanup = onCleanup(@() FFTWBackend.restoreMatlabState(state));
            try
                capabilities = FFTWBackend.inspectCapabilities(context);
            catch exception
                capabilities = FFTWBackend.emptyCapabilities(context);
                capabilities.reason = FFTWBackend.reason("inspection-failed","inspection","",string(exception.identifier),string(exception.message));
                capabilities.status = "unavailable";
            end
            clear stateCleanup
        end

        function capabilities = buildWithContext(context)
            state = FFTWBackend.captureMatlabState;
            stateCleanup = onCleanup(@() FFTWBackend.restoreMatlabState(state));
            buildStart = tic;
            startedAtUTC = FFTWBackend.utcTimestamp;
            if FFTWBackend.hasLivePlans("fftw_r2c") || FFTWBackend.hasLivePlans("fftw_r2r")
                capabilities = FFTWBackend.emptyCapabilities(context);
                capabilities.reason = FFTWBackend.reason("module-in-use","build","","FFTWBackend:ModuleInUse","An FFTW MEX module is locked by a live transform. Delete all transforms before rebuilding.");
                capabilities.build.attempted = true;
                capabilities.build.succeeded = false;
                capabilities.build.startedAtUTC = startedAtUTC;
                capabilities.build.completedAtUTC = FFTWBackend.utcTimestamp;
                capabilities.build.durationSeconds = toc(buildStart);
                capabilities.build.reason = capabilities.reason;
                clear stateCleanup
                return
            end
            capabilities = FFTWBackend.inspectCapabilities(context);
            capabilities.build.attempted = true;
            capabilities.build.startedAtUTC = startedAtUTC;
            installStarted = false;
            destinations = FFTWBackend.modulePaths(context.sourceDirectory,context.mexExtension);
            backups = strings(1,2);
            hadOriginal = false(1,2);
            try
                FFTWBackend.requireBuildPreflight(capabilities,context);
                if context.failureInjection == "compile-failed"
                    error('FFTWBackend:InjectedCompileFailure','Injected compilation failure.');
                end

                stageDirectory = string(tempname);
                [created,message] = mkdir(stageDirectory);
                if ~created
                    error('FFTWBackend:StagingFailed','Unable to create the staging directory: %s',message);
                end
                stageCleanup = onCleanup(@() FFTWBackend.removeTemporaryDirectory(stageDirectory));
                FFTWBackend.compileModules(context,stageDirectory);
                FFTWBackend.clearMexModules;

                stagedContext = context;
                stagedContext.moduleDirectory = stageDirectory;
                stagedCapabilities = FFTWBackend.inspectCapabilities(stagedContext);
                if ~stagedCapabilities.isComplete
                    error('FFTWBackend:StagedValidationFailed','The staged FFTW modules failed identity, alignment, or numerical validation: %s',stagedCapabilities.reason.message);
                end
                FFTWBackend.clearMexModules;

                backupDirectory = fullfile(stageDirectory,"backup");
                mkdir(backupDirectory);
                for iModule = 1:2
                    hadOriginal(iModule) = isfile(destinations(iModule));
                    backups(iModule) = fullfile(backupDirectory,"module-"+iModule+"."+context.mexExtension);
                    if hadOriginal(iModule)
                        [copied,message] = copyfile(destinations(iModule),backups(iModule),'f');
                        if ~copied, error('FFTWBackend:BackupFailed','Unable to back up %s: %s',destinations(iModule),message); end
                    end
                end

                installStarted = true;
                stagedPaths = FFTWBackend.modulePaths(stageDirectory,context.mexExtension);
                [copied,message] = copyfile(stagedPaths(1),destinations(1),'f');
                if ~copied, error('FFTWBackend:InstallFailed','Unable to install fftw_r2c: %s',message); end
                if context.failureInjection == "install-failed"
                    error('FFTWBackend:InjectedInstallFailure','Injected installation failure.');
                end
                [copied,message] = copyfile(stagedPaths(2),destinations(2),'f');
                if ~copied, error('FFTWBackend:InstallFailed','Unable to install fftw_r2r: %s',message); end
                FFTWBackend.clearMexModules;
                rehash

                finalContext = context;
                finalContext.moduleDirectory = context.sourceDirectory;
                finalCapabilities = FFTWBackend.inspectCapabilities(finalContext);
                if context.failureInjection == "final-validation-failed"
                    error('FFTWBackend:InjectedFinalValidationFailure','Injected final validation failure.');
                end
                if ~finalCapabilities.isComplete
                    error('FFTWBackend:FinalValidationFailed','The installed FFTW modules failed final validation: %s',finalCapabilities.reason.message);
                end
                capabilities = finalCapabilities;
                capabilities.build.attempted = true;
                capabilities.build.succeeded = true;
                capabilities.build.installed = true;
                capabilities.build.startedAtUTC = startedAtUTC;
                capabilities.build.reason = FFTWBackend.emptyReason;
                clear stageCleanup
            catch exception
                if installStarted
                    FFTWBackend.clearMexModules;
                    FFTWBackend.restoreInstalledModules(destinations,backups,hadOriginal);
                    rehash
                end
                capabilities = FFTWBackend.inspectCapabilities(context);
                capabilities.build.attempted = true;
                capabilities.build.succeeded = false;
                capabilities.build.installed = false;
                capabilities.build.startedAtUTC = startedAtUTC;
                capabilities.build.reason = FFTWBackend.exceptionReason(exception,"build");
            end
            capabilities.build.durationSeconds = toc(buildStart);
            capabilities.build.completedAtUTC = FFTWBackend.utcTimestamp;
            clear stateCleanup
        end
    end

    methods (Static, Access=private)
        function provider = bundledProvider(context)
            provider.id = "matlab-bundled";
            provider.type = "bundled";
            provider.origin = "active MATLAB installation";
            provider.isDefault = true;
            provider.distributionPolicy = "build locally; distribute no MEX or FFTW binary";
            provider.supportedReleases = "2026a";
            provider.supportedArchitectures = "maca64";
            provider.libraryPath = fullfile(context.matlabRoot,"bin",context.architecture,"libmwfftw3.3.dylib");
            provider.mexArguments = "-R2018a";
            provider.linkArguments = provider.libraryPath;
        end

        function capabilities = inspectCapabilities(context)
            if FFTWBackend.canonicalPath(context.moduleDirectory) ~= FFTWBackend.canonicalPath(context.sourceDirectory)
                previousPath = path;
                pathCleanup = onCleanup(@() path(previousPath));
                addpath(context.moduleDirectory,'-begin');
            end
            capabilities = FFTWBackend.emptyCapabilities(context);

            if context.release ~= context.provider.supportedReleases
                capabilities.reason = FFTWBackend.reason("unsupported-release","platform","","",sprintf('MATLAB release %s is unsupported; this backend currently supports R2026a.',context.release));
                capabilities.build.isPossible = false;
                return
            end
            if context.architecture ~= context.provider.supportedArchitectures
                capabilities.reason = FFTWBackend.reason("unsupported-architecture","platform","","",sprintf('Architecture %s is unsupported; this backend currently supports maca64.',context.architecture));
                capabilities.build.isPossible = false;
                return
            end
            capabilities.platform.isSupported = true;
            if context.failureInjection == "compiler-unavailable"
                capabilities.compiler.isAvailable = false;
            end
            if ~capabilities.compiler.isAvailable
                capabilities.reason = FFTWBackend.reason("compiler-unavailable","compiler","","","No selected MATLAB C++ MEX compiler is available.");
                capabilities.build.isPossible = false;
            end

            expectedLibrary = FFTWBackend.canonicalPath(context.provider.libraryPath);
            capabilities.library.expectedPath = expectedLibrary;
            libraryExists = isfile(expectedLibrary) && context.failureInjection ~= "library-missing";
            capabilities.library.exists = libraryExists;
            if ~libraryExists
                capabilities.reason = FFTWBackend.reason("bundled-library-missing","library","","",sprintf('The active MATLAB FFTW library was not found at %s.',expectedLibrary));
                capabilities.build.isPossible = false;
                return
            end

            capabilities.modules.r2c = FFTWBackend.inspectModule("fftw_r2c",context,expectedLibrary);
            capabilities.modules.r2r = FFTWBackend.inspectModule("fftw_r2r",context,expectedLibrary);
            if capabilities.modules.r2c.isAvailable
                [r2c,c2r] = FFTWBackend.testR2CModule(context);
                capabilities.features.r2c = r2c;
                capabilities.features.c2r = c2r;
            else
                capabilities.features.r2c.reason = capabilities.modules.r2c.reason;
                capabilities.features.c2r.reason = capabilities.modules.r2c.reason;
            end
            if capabilities.modules.r2r.isAvailable
                [dct1,dst1] = FFTWBackend.testR2RModule(context);
                capabilities.features.dct1 = dct1;
                capabilities.features.dst1 = dst1;
            else
                capabilities.features.dct1.reason = capabilities.modules.r2r.reason;
                capabilities.features.dst1.reason = capabilities.modules.r2r.reason;
            end

            flags = [capabilities.features.r2c.isAvailable capabilities.features.c2r.isAvailable capabilities.features.dct1.isAvailable capabilities.features.dst1.isAvailable];
            capabilities.isAvailable = any(flags);
            capabilities.isComplete = all(flags);
            if capabilities.isComplete
                capabilities.status = "available";
                capabilities.reason = FFTWBackend.emptyReason;
            elseif capabilities.isAvailable
                capabilities.status = "partial";
                capabilities.reason = FFTWBackend.firstFeatureReason(capabilities.features);
            else
                capabilities.status = "unavailable";
                capabilities.reason = FFTWBackend.firstFeatureReason(capabilities.features);
            end
            capabilities.library.resolvedPath = FFTWBackend.commonLibraryPath(capabilities.modules);
            capabilities.library.version = FFTWBackend.commonLibraryVersion(capabilities.modules);
            capabilities.library.identityValidated = capabilities.modules.r2c.identityValidated && capabilities.modules.r2r.identityValidated;
            capabilities.provider.identityValidated = capabilities.library.identityValidated;
            capabilities.build.isRequired = ~capabilities.isComplete;
            capabilities.build.isPossible = capabilities.platform.isSupported && capabilities.compiler.isAvailable && capabilities.library.exists;
            clear pathCleanup
        end

        function module = inspectModule(name,context,expectedLibrary)
            module = FFTWBackend.emptyModule(name);
            expectedPath = fullfile(context.moduleDirectory,name+"."+context.mexExtension);
            module.expectedPath = FFTWBackend.canonicalPath(expectedPath);
            missingInjection = context.failureInjection == name+"-missing";
            if ~isfile(expectedPath) || missingInjection
                module.reason = FFTWBackend.reason("mex-missing","module",name,"",sprintf('%s is not built at %s.',name,expectedPath));
                return
            end
            resolvedModule = string(which(name));
            if context.failureInjection == "unexpected-module-path"
                resolvedModule = fullfile(tempdir,name+"."+context.mexExtension);
            end
            module.path = FFTWBackend.canonicalPath(resolvedModule);
            module.exists = true;
            if module.path ~= module.expectedPath
                module.reason = FFTWBackend.reason("unexpected-module-path","module",name,"",sprintf('%s resolved to %s instead of %s.',name,module.path,module.expectedPath));
                return
            end
            try
                if context.failureInjection == "missing-symbol"
                    error('FFTWBackend:MissingSymbol','Injected missing FFTW symbol.');
                end
                [versionText,libraryPath] = feval(name,'info');
                module.loaded = true;
                module.libraryVersion = string(versionText);
                module.libraryPath = FFTWBackend.canonicalPath(string(libraryPath));
                if context.failureInjection == "library-mismatch"
                    module.libraryPath = fullfile(tempdir,"libfftw3.dylib");
                end
                if module.libraryPath ~= expectedLibrary
                    module.reason = FFTWBackend.reason("library-mismatch","identity",name,"",sprintf('%s loaded %s instead of the active MATLAB library %s.',name,module.libraryPath,expectedLibrary));
                    return
                end
                module.identityValidated = true;
                module.isAvailable = true;
                module.reason = FFTWBackend.emptyReason;
            catch exception
                code = "mex-load-failed";
                if string(exception.identifier) == "FFTWBackend:MissingSymbol", code = "missing-symbol"; end
                module.reason = FFTWBackend.reason(code,"load",name,string(exception.identifier),string(exception.message));
            end
        end

        function [r2c,c2r] = testR2CModule(context)
            r2c = FFTWBackend.emptyFeature("r2c","fftw_r2c");
            c2r = FFTWBackend.emptyFeature("c2r","fftw_r2c");
            try
                input = reshape(sin((1:(8*6*3))/7),[8 6 3]);
                transform = RealToComplexTransform(size(input),dims=[2 1],planner="estimate",nCores=1,alignmentMode="unaligned",plannerTimeLimitSeconds=1);
                transformCleanup = onCleanup(@() delete(transform));
                spectrum = transform.transformForward(input);
                expected = fft(fft(input,[],2),[],1);
                expected = expected(1:5,:,:);
                r2c.maximumRelativeError = FFTWBackend.relativeError(spectrum,expected);
                inverse = transform.scaleFactor*transform.transformBack(spectrum);
                c2r.maximumRelativeError = FFTWBackend.relativeError(inverse,input);
                clear transformCleanup

                matched = RealToComplexTransform(size(input),dims=[2 1],planner="estimate",nCores=1,alignmentMode="matched",plannerTimeLimitSeconds=1);
                matchedCleanup = onCleanup(@() delete(matched));
                matchedSpectrum = matched.transformForward(input);
                matchedError = FFTWBackend.relativeError(matchedSpectrum,expected);
                [matchedAccepted,mismatchRejected,unalignedAccepted] = fftw_r2c('alignmentSelfTest');
                alignmentPassed = matchedAccepted && mismatchRejected && unalignedAccepted && matchedError <= context.errorTolerance;
                if context.failureInjection == "alignment-failed", alignmentPassed = false; end
                r2c.alignmentPassed = alignmentPassed;
                c2r.alignmentPassed = alignmentPassed;
                clear matchedCleanup

                if context.failureInjection == "r2c-self-test-failed", r2c.maximumRelativeError = Inf; end
                if context.failureInjection == "c2r-self-test-failed", c2r.maximumRelativeError = Inf; end
                r2c = FFTWBackend.finalizeFeature(r2c,context.errorTolerance);
                c2r = FFTWBackend.finalizeFeature(c2r,context.errorTolerance);
            catch exception
                reason = FFTWBackend.reason("self-test-failed","self-test","fftw_r2c",string(exception.identifier),string(exception.message));
                r2c.reason = reason;
                c2r.reason = reason;
            end
        end

        function [dct1,dst1] = testR2RModule(context)
            dct1 = FFTWBackend.emptyFeature("dct1","fftw_r2r");
            dst1 = FFTWBackend.emptyFeature("dst1","fftw_r2r");
            try
                n = 9;
                values = reshape(sin((1:(n*3))/5),[n 3])+1i*reshape(cos((1:(n*3))/11),[n 3]);
                dct = RealToRealTransform(size(values),dims=1,transform="cosine",dataType="complex",planner="estimate",nCores=1,alignmentMode="unaligned",plannerTimeLimitSeconds=1);
                dctCleanup = onCleanup(@() delete(dct));
                dctCoefficients = dct.transformForward(values);
                [dctForwardMatrix,dctBackMatrix] = FFTWBackend.referenceMatrices(n,"cosine");
                dctExpected = dctForwardMatrix*values;
                dctForwardError = FFTWBackend.relativeError(dctCoefficients,dctExpected);
                dctBack = dct.transformBack(dctCoefficients);
                dctBackExpected = dctBackMatrix*dctCoefficients;
                dct1.maximumRelativeError = max(dctForwardError,FFTWBackend.relativeError(dctBack,dctBackExpected));
                clear dctCleanup

                sineValues = values;
                sineValues(1,:) = 3+2i;
                sineValues(end,:) = -4+5i;
                dst = RealToRealTransform(size(values),dims=1,transform="sine",dataType="complex",planner="estimate",nCores=1,alignmentMode="unaligned",plannerTimeLimitSeconds=1);
                dstCleanup = onCleanup(@() delete(dst));
                dstCoefficients = dst.transformForward(sineValues);
                [dstForwardMatrix,dstBackMatrix] = FFTWBackend.referenceMatrices(n,"sine");
                dstExpected = dstForwardMatrix*sineValues;
                dstForwardError = FFTWBackend.relativeError(dstCoefficients,dstExpected);
                dstBack = dst.transformBack(dstCoefficients);
                dstBackExpected = dstBackMatrix*dstCoefficients;
                endpointPassed = all(dstBack([1 end],:) == 0,'all');
                dst1.maximumRelativeError = max(dstForwardError,FFTWBackend.relativeError(dstBack,dstBackExpected));
                clear dstCleanup

                matched = RealToRealTransform([n 3],dims=1,transform="cosine",planner="estimate",nCores=1,alignmentMode="matched",plannerTimeLimitSeconds=1);
                matchedCleanup = onCleanup(@() delete(matched));
                matchedError = FFTWBackend.relativeError(matched.transformForward(real(values)),dctForwardMatrix*real(values));
                [matchedAccepted,mismatchRejected,unalignedAccepted] = fftw_r2r('alignmentSelfTest');
                alignmentPassed = matchedAccepted && mismatchRejected && unalignedAccepted && matchedError <= context.errorTolerance;
                if context.failureInjection == "alignment-failed", alignmentPassed = false; end
                dct1.alignmentPassed = alignmentPassed;
                dst1.alignmentPassed = alignmentPassed && endpointPassed;
                clear matchedCleanup

                if context.failureInjection == "dct1-self-test-failed", dct1.maximumRelativeError = Inf; end
                if context.failureInjection == "dst1-self-test-failed", dst1.maximumRelativeError = Inf; end
                dct1 = FFTWBackend.finalizeFeature(dct1,context.errorTolerance);
                dst1 = FFTWBackend.finalizeFeature(dst1,context.errorTolerance);
            catch exception
                reason = FFTWBackend.reason("self-test-failed","self-test","fftw_r2r",string(exception.identifier),string(exception.message));
                dct1.reason = reason;
                dst1.reason = reason;
            end
        end

        function feature = finalizeFeature(feature,tolerance)
            feature.selfTestPassed = isfinite(feature.maximumRelativeError) && feature.maximumRelativeError <= tolerance;
            feature.isAvailable = feature.selfTestPassed && feature.alignmentPassed;
            if feature.isAvailable
                feature.reason = FFTWBackend.emptyReason;
            elseif ~feature.alignmentPassed
                feature.reason = FFTWBackend.reason("alignment-failed","alignment",feature.module,"","FFTW alignment validation failed.");
            else
                feature.reason = FFTWBackend.reason("self-test-failed","self-test",feature.module,"",sprintf('%s relative error %.3g exceeds %.3g.',feature.id,feature.maximumRelativeError,tolerance));
            end
        end

        function compileModules(context,stageDirectory)
            includeArgument = "-I"+context.sourceDirectory;
            mex(context.provider.mexArguments,'-outdir',stageDirectory,'-output','fftw_r2c',includeArgument,fullfile(context.sourceDirectory,'fftw_r2c.cpp'),context.provider.linkArguments);
            mex(context.provider.mexArguments,'-outdir',stageDirectory,'-output','fftw_r2r',includeArgument,fullfile(context.sourceDirectory,'fftw_r2r.cpp'),context.provider.linkArguments);
            rehash
        end

        function requireBuildPreflight(capabilities,context)
            if ~capabilities.platform.isSupported
                error('FFTWBackend:UnsupportedPlatform','%s',capabilities.reason.message);
            end
            if ~capabilities.compiler.isAvailable
                error('FFTWBackend:CompilerUnavailable','No selected MATLAB C++ MEX compiler is available.');
            end
            if ~isfile(context.provider.libraryPath) || context.failureInjection == "library-missing"
                error('FFTWBackend:LibraryMissing','The active MATLAB FFTW library is unavailable.');
            end
        end

        function restoreInstalledModules(destinations,backups,hadOriginal)
            for iModule = 1:2
                if hadOriginal(iModule)
                    copyfile(backups(iModule),destinations(iModule),'f');
                elseif isfile(destinations(iModule))
                    delete(destinations(iModule));
                end
            end
        end

        function removeTemporaryDirectory(path)
            if isfolder(path), rmdir(path,'s'); end
        end

        function clearMexModules()
            clear fftw_r2c fftw_r2r
            rehash
        end

        function value = hasLivePlans(name)
            value = false;
            if exist(name,'file') ~= 3, return, end
            try
                lifetime = feval(name,'lifetime');
                value = lifetime(3) > 0;
            catch
                value = mislocked(name);
            end
        end

        function paths = modulePaths(directory,extension)
            paths = [fullfile(directory,"fftw_r2c."+extension) fullfile(directory,"fftw_r2r."+extension)];
        end

        function capabilities = emptyCapabilities(context)
            capabilities.schemaVersion = "1.0.0";
            capabilities.status = "unavailable";
            capabilities.isAvailable = false;
            capabilities.isComplete = false;
            capabilities.platform.matlabVersion = context.matlabVersion;
            capabilities.platform.matlabRelease = context.release;
            capabilities.platform.architecture = context.architecture;
            capabilities.platform.mexExtension = context.mexExtension;
            capabilities.platform.isSupported = false;
            capabilities.provider.id = context.provider.id;
            capabilities.provider.type = context.provider.type;
            capabilities.provider.origin = context.provider.origin;
            capabilities.provider.isDefault = context.provider.isDefault;
            capabilities.provider.distributionPolicy = context.provider.distributionPolicy;
            capabilities.provider.identityValidated = false;
            capabilities.compiler = FFTWBackend.compilerRecord;
            capabilities.library.expectedPath = FFTWBackend.canonicalPath(context.provider.libraryPath);
            capabilities.library.exists = isfile(context.provider.libraryPath);
            capabilities.library.resolvedPath = "";
            capabilities.library.version = "";
            capabilities.library.identityValidated = false;
            capabilities.modules.r2c = FFTWBackend.emptyModule("fftw_r2c");
            capabilities.modules.r2r = FFTWBackend.emptyModule("fftw_r2r");
            capabilities.features.r2c = FFTWBackend.emptyFeature("r2c","fftw_r2c");
            capabilities.features.c2r = FFTWBackend.emptyFeature("c2r","fftw_r2c");
            capabilities.features.dct1 = FFTWBackend.emptyFeature("dct1","fftw_r2r");
            capabilities.features.dst1 = FFTWBackend.emptyFeature("dst1","fftw_r2r");
            capabilities.eligibility = FFTWBackend.eligibility;
            capabilities.build.isRequired = true;
            capabilities.build.isPossible = capabilities.compiler.isAvailable && capabilities.library.exists;
            capabilities.build.attempted = false;
            capabilities.build.succeeded = false;
            capabilities.build.installed = false;
            capabilities.build.startedAtUTC = "";
            capabilities.build.completedAtUTC = "";
            capabilities.build.durationSeconds = NaN;
            capabilities.build.reason = FFTWBackend.emptyReason;
            capabilities.reason = FFTWBackend.emptyReason;
        end

        function compiler = compilerRecord()
            compiler.name = "";
            compiler.version = "";
            compiler.isAvailable = false;
            try
                configuration = mex.getCompilerConfigurations('C++','Selected');
                if ~isempty(configuration)
                    compiler.name = string(configuration(1).Name);
                    compiler.version = string(configuration(1).Version);
                    compiler.isAvailable = true;
                end
            catch
            end
        end

        function module = emptyModule(name)
            module.name = name;
            module.expectedPath = "";
            module.path = "";
            module.exists = false;
            module.loaded = false;
            module.libraryPath = "";
            module.libraryVersion = "";
            module.identityValidated = false;
            module.isAvailable = false;
            module.reason = FFTWBackend.emptyReason;
        end

        function feature = emptyFeature(id,module)
            feature.id = id;
            feature.module = module;
            feature.isAvailable = false;
            feature.selfTestPassed = false;
            feature.alignmentPassed = false;
            feature.maximumRelativeError = NaN;
            feature.reason = FFTWBackend.emptyReason;
        end

        function eligibility = eligibility()
            eligibility.horizontal.schemaVersion = "issue41-v1";
            eligibility.horizontal.isReady = true;
            eligibility.horizontal.layout = "half-x";
            eligibility.horizontal.transformDimensions = [2 1];
            eligibility.horizontal.strategy = "guru-rank2";
            eligibility.horizontal.ownership = "MATLAB-managed zero-copy forward and uniquely owned destructive inverse";
            eligibility.horizontal.gateSizes = [256 256 64; 1024 1024 30];
            eligibility.horizontal.thresholds = struct('rawForwardSpeedup',1.25,'completeForwardSpeedup',1.10,'destructiveInverseSpeedRatio',0.95,'maximumRelativeError',1e-12);
            eligibility.horizontal.sourceIssue = 41;
            eligibility.horizontal.sourceCommit = "a7aef07";
            eligibility.horizontal.sourceArtifact = "Matlab/Spectral/FFTW/results/issue41/20260808T153906007Z-maca64-r2026a/bundled-native-comparison.json";

            eligibility.realToReal.schemaVersion = "issue43-v1";
            eligibility.realToReal.policy = "exact Nz and inclusive bounded batch intervals only; no extrapolation";
            eligibility.realToReal.testedNz = [33 65 129 257 513];
            eligibility.realToReal.testedBatchCounts = [1 8320 33024 131584];
            eligibility.realToReal.speedThreshold = 1.10;
            eligibility.realToReal.maximumRelativeError = 1e-12;
            eligibility.realToReal.sourceIssue = 43;
            eligibility.realToReal.sourceCommit = "18a37ce";
            eligibility.realToReal.sourceArtifact = "Matlab/Spectral/FFTW/results/issue43/20260808T163913752Z-maca64-r2026a/real-to-real-benchmark.json";
            eligibility.realToReal.records = FFTWBackend.realToRealEligibilityRecords;
        end

        function records = realToRealEligibilityRecords()
            nzValues = [33 65 129 257 513];
            keys = {
                "real","cosine","forward"
                "real","cosine","inverse"
                "real","sine","forward"
                "real","sine","inverse"
                "complex","cosine","forward"
                "complex","cosine","inverse"
                "complex","sine","forward"
                "complex","sine","inverse"
                };
            minimums = [
                NaN NaN 131584 NaN 33024 33024 8320 131584
                NaN NaN 33024 NaN 33024 33024 8320 33024
                NaN NaN 8320 NaN 8320 8320 8320 8320
                8320 NaN 8320 8320 8320 8320 8320 8320
                8320 8320 8320 8320 8320 8320 8320 8320
                ];
            emptyIntervals = repmat(struct('minimumBatchCount',0,'maximumBatchCount',0),0,1);
            template = struct('Nz',0,'dataType',"",'transformType',"",'direction',"",'eligible',false,'intervals',emptyIntervals,'testedBatchCounts',[1 8320 33024 131584]);
            records = repmat(template,numel(nzValues)*size(keys,1),1);
            iRecord = 0;
            for iNz = 1:numel(nzValues)
                for iKey = 1:size(keys,1)
                    iRecord = iRecord+1;
                    records(iRecord).Nz = nzValues(iNz);
                    records(iRecord).dataType = keys{iKey,1};
                    records(iRecord).transformType = keys{iKey,2};
                    records(iRecord).direction = keys{iKey,3};
                    if isfinite(minimums(iNz,iKey))
                        records(iRecord).eligible = true;
                        records(iRecord).intervals = struct('minimumBatchCount',minimums(iNz,iKey),'maximumBatchCount',131584);
                    end
                end
            end
        end

        function reason = firstFeatureReason(features)
            names = ["r2c","c2r","dct1","dst1"];
            reason = FFTWBackend.emptyReason;
            for name = names
                candidate = features.(name).reason;
                if strlength(candidate.code) > 0
                    reason = candidate;
                    return
                end
            end
        end

        function path = commonLibraryPath(modules)
            paths = [modules.r2c.libraryPath modules.r2r.libraryPath];
            paths = paths(strlength(paths) > 0);
            if isempty(paths), path = ""; elseif all(paths == paths(1)), path = paths(1); else, path = "multiple"; end
        end

        function value = commonLibraryVersion(modules)
            versions = [modules.r2c.libraryVersion modules.r2r.libraryVersion];
            versions = versions(strlength(versions) > 0);
            if isempty(versions), value = ""; elseif all(versions == versions(1)), value = versions(1); else, value = "multiple"; end
        end

        function value = relativeError(actual,expected)
            value = norm(actual(:)-expected(:),inf)/max(norm(expected(:),inf),eps);
        end

        function [forward,back] = referenceMatrices(n,transformType)
            indices = 0:(n-1);
            if transformType == "cosine"
                forward = (2/(n-1))*cos(pi*(indices.')*indices/(n-1));
                forward([1 end],:) = forward([1 end],:)/2;
                forward(:,[1 end]) = forward(:,[1 end])/2;
                forward(1,:) = 2*forward(1,:);
                back = cos(pi*(indices.')*indices/(n-1));
                back(:,1) = back(:,1)/2;
            else
                interior = 1:(n-2);
                forward = (2/(n-1))*sin(pi*(interior.')*indices/(n-1));
                back = sin(pi*(indices.')*interior/(n-1));
            end
        end

        function path = canonicalPath(path)
            if strlength(string(path)) == 0, path = ""; return, end
            try
                path = string(java.io.File(char(path)).getCanonicalPath());
            catch
                path = string(path);
            end
        end

        function reason = exceptionReason(exception,stage)
            switch string(exception.identifier)
                case "FFTWBackend:UnsupportedPlatform"
                    code = "unsupported-platform";
                case "FFTWBackend:CompilerUnavailable"
                    code = "compiler-unavailable";
                case "FFTWBackend:LibraryMissing"
                    code = "bundled-library-missing";
                case "FFTWBackend:InjectedCompileFailure"
                    code = "compile-failed";
                case "FFTWBackend:StagingFailed"
                    code = "staging-failed";
                case "FFTWBackend:StagedValidationFailed"
                    code = "staged-validation-failed";
                case "FFTWBackend:BackupFailed"
                    code = "backup-failed";
                case {"FFTWBackend:InstallFailed","FFTWBackend:InjectedInstallFailure"}
                    code = "install-failed";
                case {"FFTWBackend:FinalValidationFailed","FFTWBackend:InjectedFinalValidationFailure"}
                    code = "final-validation-failed";
                otherwise
                    code = stage+"-failed";
            end
            reason = FFTWBackend.reason(code,stage,"",string(exception.identifier),string(exception.message));
        end

        function value = reason(code,stage,module,identifier,message)
            value.code = string(code);
            value.stage = string(stage);
            value.module = string(module);
            value.identifier = string(identifier);
            value.message = string(message);
        end

        function value = emptyReason()
            value = FFTWBackend.reason("","","","","");
        end

        function value = utcTimestamp()
            value = string(datetime('now',TimeZone='UTC',Format="yyyy-MM-dd'T'HH:mm:ss.SSS'Z'"));
        end

        function state = captureMatlabState()
            state.directory = pwd;
            state.path = path;
            state.planner = string(fftw('planner'));
            state.wisdom = fftw('dwisdom');
            state.threads = maxNumCompThreads;
            state.random = rng;
        end

        function restoreMatlabState(state)
            try
                cd(state.directory);
            catch
            end
            try
                path(state.path);
            catch
            end
            try
                fftw('dwisdom',[]);
                fftw('dwisdom',state.wisdom);
            catch
            end
            try
                fftw('planner',char(state.planner));
            catch
            end
            try
                maxNumCompThreads(state.threads);
            catch
            end
            try
                rng(state.random);
            catch
            end
        end
    end
end
