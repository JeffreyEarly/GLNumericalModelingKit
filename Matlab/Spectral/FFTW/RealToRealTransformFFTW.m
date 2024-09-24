classdef RealToRealTransformFFTW < handle
    properties
        plan
        scaleFactor = 1;
        outp
    end

    methods
        function self = RealToRealTransformFFTW(sz,options)
            arguments
                sz
                options.dim double = 1
                options.transform char {mustBeMember(options.transform,["cosine","sine"])}
                options.planner char {mustBeMember(options.planner,["estimate","measure","patient","exhaustive"])} = "estimate"
            end
            if ~isfield(options,'transform')
                error('You must specify a transform type.');
            end

            if not(libisloaded('libmwfftw3'))
                addpath(fullfile(matlabroot,'bin','maca64'))
                addpath(fullfile(matlabroot,'bin','maci64'))
                loadlibrary('libmwfftw3','fftw3.h')
            end

            % The logic here is complicated by the fact that we might have
            % singleton dimension, which we need to ignore.
            lastDim = 0;
            iDim = 0;
            for iSize=1:length(sz)
                if sz(iSize)==1
                    continue;
                end
                iDim = iDim+1;
                if iDim==1
                    alldims(iDim).n = sz(iSize);
                    alldims(iDim).is = 1;
                    alldims(iDim).os = 1;
                else
                    alldims(iDim).n = sz(iSize);
                    alldims(iDim).is = alldims(lastDim).is*alldims(lastDim).n;
                    alldims(iDim).os = alldims(lastDim).os*alldims(lastDim).n;
                end
                lastDim = iDim;
            end

            rank = 0;
            howmany_rank = 0;
            for iDim=1:length(alldims)
                if iDim==options.dim
                    rank = rank + 1;
                    dims(rank).n = alldims(iDim).n;
                    dims(rank).is = alldims(iDim).is;
                    dims(rank).os = alldims(iDim).os;
                    self.scaleFactor = 1/(dims(rank).n -1);
                else
                    howmany_rank = howmany_rank + 1;
                    howmany_dims(howmany_rank).n = alldims(iDim).n;
                    howmany_dims(howmany_rank).is = alldims(iDim).is;
                    howmany_dims(howmany_rank).os = alldims(iDim).os;
                end
            end

% See GLBasisTransformationOperation for how to actually do this.

% #define FFTW_MEASURE (0U)
% #define FFTW_DESTROY_INPUT (1U << 0)
% #define FFTW_UNALIGNED (1U << 1)
% #define FFTW_CONSERVE_MEMORY (1U << 2)
% #define FFTW_EXHAUSTIVE (1U << 3) /* NO_EXHAUSTIVE is default */
% #define FFTW_PRESERVE_INPUT (1U << 4) /* cancels FFTW_DESTROY_INPUT */
% #define FFTW_PATIENT (1U << 5) /* IMPATIENT is default */
% #define FFTW_ESTIMATE (1U << 6)
% #define FFTW_WISDOM_ONLY (1U << 21)

% FFTW_R2HC=0, FFTW_HC2R=1, FFTW_DHT=2,
% FFTW_REDFT00=3, FFTW_REDFT01=4, FFTW_REDFT10=5, FFTW_REDFT11=6,
% FFTW_RODFT00=7, FFTW_RODFT01=8, FFTW_RODFT10=9, FFTW_RODFT11=10

            switch options.transform
                case "cosine"
                    transformKind = 3;
                case "sine"
                    transformKind = 7;
            end
            switch options.planner
                case "estimate"
                    planner = bitshift(1,6);
                case "measure"
                    planner = 0;
                case "patient"
                    planner = bitshift(1,5);
                case "exhaustive"
                    planner = bitshift(1,3);
            end

            in = zeros(sz);
            out = zeros(sz);
            self.plan = calllib('libmwfftw3','fftw_plan_guru_r2r',rank,dims,howmany_rank,howmany_dims,in,out,transformKind,planner);
            self.outp = libpointer('doublePtr',zeros(sz));
        end

        function f = transformBack(self,fbar)
            
            calllib('libmwfftw3','fftw_execute_r2r', self.plan, fbar, self.outp );
            f = self.outp.Value/2;
        end

        function fbar = transformForward(self,f)
            calllib('libmwfftw3','fftw_execute_r2r', self.plan, f, self.outp );
            fbar = self.scaleFactor*self.outp.Value;
        end
    end
end