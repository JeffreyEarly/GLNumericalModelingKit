classdef RealToRealTransformMexFFTW < RealToRealTransform
    % Compatibility name for RealToRealTransform.
    %
    % New code should construct RealToRealTransform directly. This class
    % translates the legacy dim option and delegates all behavior to the
    % production bundled-FFTW implementation.

    methods
        function self = RealToRealTransformMexFFTW(sz,options)
            arguments
                sz (1,:) double {mustBeInteger,mustBePositive}
                options.dim (1,1) double {mustBeInteger,mustBePositive} = 1
                options.transform (1,1) string {mustBeMember(options.transform,["cosine","sine"])}
                options.dataType (1,1) string {mustBeMember(options.dataType,["real","complex"])} = "real"
                options.planner (1,1) string {mustBeMember(options.planner,["estimate","measure","patient","exhaustive"])} = "measure"
                options.nCores (1,1) double {mustBeInteger,mustBePositive} = maxNumCompThreads
                options.alignmentMode (1,1) string {mustBeMember(options.alignmentMode,["matched","unaligned"])} = "unaligned"
                options.plannerTimeLimitSeconds (1,1) double {mustBePositive,mustBeFinite} = 10
            end
            self@RealToRealTransform(sz,dims=options.dim,transform=options.transform,dataType=options.dataType,planner=options.planner,nCores=options.nCores,alignmentMode=options.alignmentMode,plannerTimeLimitSeconds=options.plannerTimeLimitSeconds);
        end
    end

    methods (Static)
        function makeMexFiles(fftwlibpath)
            arguments
                fftwlibpath (1,1) string = fullfile(matlabroot,'bin',computer('arch'),'libmwfftw3.3.dylib')
            end
            RealToRealTransform.makeMexFiles(fftwlibpath);
        end
    end
end
