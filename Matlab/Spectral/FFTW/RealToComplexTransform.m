classdef RealToComplexTransform < handle
    properties
        realSize
        complexSize
        scaleFactor = 1;
        plan
    end

    % properties (Access=private)
    %     plan
    % end

    methods
        function self = RealToComplexTransform(sz,options)
            arguments
                sz
                options.dims double = 1
                options.planner char {mustBeMember(options.planner,["estimate","measure","patient","exhaustive"])} = "estimate"
                options.nCores = 1
            end
            self.realSize = sz;

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
                    planner = bitshift(1,6);
                case "measure"
                    planner = 0;
                case "patient"
                    planner = bitshift(1,5);
                case "exhaustive"
                    planner = bitshift(1,3);
            end

            % self.plan = fftw_plan_guru_dft_r2c(sz, options.dims, options.nCores, self.planner);
            [self.plan, self.complexSize, self.scaleFactor] = fftw_plan_guru_dft_r2c(sz, options.dims, options.nCores, planner);
        end

        function xbar = transformForward(self,x)
            xbar = fftw_execute_dft_r2c(self.plan,x);
        end

        function x = transformBack(self,xbar)
            x = fftw_execute_dft_c2r(self.plan,xbar);
        end
    end

    methods (Static)
        function makeMexFiles(fftwlibpath)
            arguments
                fftwlibpath = fullfile(matlabroot,'bin',computer('arch'),'libmwfftw3.3.dylib')
            end
            ipath = ['-I' fullfile(matlabroot,'extern','include')];

            mex(ipath,'fftw_plan_guru_dft_r2c.cpp',fftwlibpath);
            mex(ipath,'fftw_execute_dft_r2c.cpp',fftwlibpath);
            mex(ipath,'fftw_execute_dft_c2r.cpp',fftwlibpath);
        end
    end
end