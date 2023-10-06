classdef CosineTransformFFTW < handle
    properties
        plan
    end

    methods
        function self = CosineTransformFFTW(sz,options)
            arguments
                sz
                options.dim double = 1
            end
            if not(libisloaded('libmwfftw3'))
                addpath(fullfile(matlabroot,'bin','maca64'))
                loadlibrary('libmwfftw3','fftw3.h')
            end

            % Structure indicate the dimensions to be transformed
            rank = 1;
            dims.n = sz(1);
            dims.is = 1;
            dims.os = 1;

            % Structure indicating the dimensions to be looped over
            howmany_rank = 0;
            howmany_dims.n = 0;
            howmany_dims.is = 1;
            howmany_dims.os = 1;

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

            transformKind = 3; % FFTW_REDFT00
            planner = bitshift(1,6);

            in = zeros(sz);
            out = zeros(sz);
            self.plan = calllib('libmwfftw3','fftw_plan_guru_r2r',rank,dims,howmany_rank,howmany_dims,in,out,transformKind,planner);
        end

        function f = transformBack(self,fbar)
            outp = libpointer('doublePtr',zeros(size(fbar)));
            calllib('libmwfftw3','fftw_execute_r2r', self.plan, fbar, outp );
            f = outp.Value/2;
        end

        function fbar = transformForward(self,f)
            outp = libpointer('doublePtr',zeros(size(f)));
            calllib('libmwfftw3','fftw_execute_r2r', self.plan, f, outp );
            fbar = outp.Value/(length(f)-1);
        end
    end
end