classdef RealToRealTransform < handle
    properties
        sz
        plan
        scaleFactor = 1;
        transformKind
        planner
    end

    methods
        function self = RealToRealTransform(sz,options)
            arguments
                sz
                options.dim double = 1
                options.transform char {mustBeMember(options.transform,["cosine","sine"])}
                options.planner char {mustBeMember(options.planner,["estimate","measure","patient","exhaustive"])} = "measure"
                options.nCores = 1
            end
            self.sz = sz;
            if ~isfield(options,'transform')
                error('You must specify a transform type.');
            end
            % FFTW_R2HC=0, FFTW_HC2R=1, FFTW_DHT=2,
            % FFTW_REDFT00=3, FFTW_REDFT01=4, FFTW_REDFT10=5, FFTW_REDFT11=6,
            % FFTW_RODFT00=7, FFTW_RODFT01=8, FFTW_RODFT10=9, FFTW_RODFT11=10
            switch options.transform
                case "cosine"
                    self.transformKind = 3;
                case "sine"
                    self.transformKind = 7;
            end

            % #define FFTW_MEASURE (0U)
            % #define FFTW_DESTROY_INPUT (1U << 0)
            % #define FFTW_UNALIGNED (1U << 1)
            % #define FFTW_CONSERVE_MEMORY (1U << 2)
            % #define FFTW_EXHAUSTIVE (1U << 3) /* NO_EXHAUSTIVE is default */
            % #define FFTW_PRESERVE_INPUT (1U << 4) /* cancels FFTW_DESTROY_INPUT */
            % #define FFTW_PATIENT (1U << 5) /* IMPATIENT is default */
            % #define FFTW_ESTIMATE (1U << 6)
            % #define FFTW_WISDOM_ONLY (1U << 21)
            switch options.planner
                case "estimate"
                    self.planner = bitshift(1,6);
                case "measure"
                    self.planner = 0;
                case "patient"
                    self.planner = bitshift(1,5);
                case "exhaustive"
                    self.planner = bitshift(1,3);
            end
        end
    end

    methods (Abstract)
        f = transformBack(self,fbar);
        fbar = transformForward(self,f);
    end
end