classdef RealToComplexTransform < handle
    properties
        realSize
        complexSize (1,:) double 
        scaleFactor = 1;
        scratch
    end

    properties (Access=private)
        plan
    end

    methods
        function self = RealToComplexTransform(sz,options)
            arguments
                sz
                options.dims double = 1
                options.planner char {mustBeMember(options.planner,["estimate","measure","patient","exhaustive"])} = "measure"
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

            [self.plan, self.complexSize, self.scaleFactor] = fftw_dft2('create', sz, options.dims, options.nCores, planner);
            self.scratch = complex(zeros(self.complexSize));
        end

        function xbar = transformForward(self,x)
            xbar = fftw_dft2('r2c', self.plan, x);
        end

        function fbar = transformForwardIntoArray(self,f,fbar)
            fbar = fftw_dft2('r2c', self.plan,f,fbar);
        end

        function x = transformBack(self,xbar)
            % self.scratch = reshape(xbar(1:end), size(xbar));
            % x = fftw_dft2('c2r', self.plan, self.scratch);
            x = fftw_dft2('c2r', self.plan, xbar);
        end

        function x = transformBackIntoArray(self,xbar,x)
            % self.scratch = reshape(xbar(1:end), size(xbar));
            % x = fftw_dft2('c2r_inout', self.plan, self.scratch,x);
            x = fftw_dft2('c2r_inout', self.plan, xbar,x);
        end

        function delete(self)
            fftw_dft2('free', self.plan);
        end
    end

    methods (Static)
        function makeMexFiles(fftwlibpath)
            arguments
                fftwlibpath = fullfile(matlabroot,'bin',computer('arch'),'libmwfftw3.3.dylib')
            end
            ipath = ['-I' fullfile(matlabroot,'extern','include')];
            mex(ipath,'fftw_dft2.cpp',fftwlibpath)
        end
    end
end