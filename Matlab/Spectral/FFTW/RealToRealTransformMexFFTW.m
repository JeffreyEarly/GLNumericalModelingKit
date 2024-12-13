classdef RealToRealTransformMexFFTW < RealToRealTransform
    % sine and cosine transforms
    %
    % In terms of speed, it is optimal calls in the following order:
    % 1. in-place transforms, e.g., x = dct.transformBack(x)
    % 2. preallocated out-of-place transforms, e.g.,
    %   xout = zeros(size(x));
    %   xout = dct.transformBack(x);
    % 3. out-of-place transforms, y = dct.transformBack(x);
    %
    % Option 3 is significantly slower because it must allocate memory.
    methods
        function self = RealToRealTransformMexFFTW(sz,options)
            arguments
                sz
                options.dim double = 1
                options.transform char {mustBeMember(options.transform,["cosine","sine"])}
                options.planner char {mustBeMember(options.planner,["estimate","measure","patient","exhaustive"])} = "measure"
                options.nCores = 1
            end

            scArgs = namedargs2cell(options);
            self@RealToRealTransform(sz,scArgs{:});

            % self.plan = create_dct_plan_mex(self.dims, self.howmany_dims, options.nCores, self.transformKind, self.planner);
            self.plan = fftw_plan_guru_r2r(sz, options.dim, options.nCores, self.transformKind, self.planner);
        end

        function x = transformBack(self,x)
            x = execute_dct_plan_mex(self.plan,x)/2;
        end

        function f = transformBackIntoArray(self,fbar,f)     
            f = execute_dct_plan_inout_mex(self.plan,fbar,f)/2;
        end

        function x = transformForward(self,x)
            x = self.scaleFactor*execute_dct_plan_mex(self.plan,x);
        end

        function fbar = transformForwardIntoArray(self,f,fbar)
            fbar = self.scaleFactor*execute_dct_plan_inout_mex(self.plan,f,fbar);
        end
    end

    methods (Static)
        function makeMexFiles(fftwlibpath)
            arguments
                fftwlibpath = fullfile(matlabroot,'bin',computer('arch'),'libmwfftw3.3.dylib')
            end
            % fftwlib = fullfile(matlabroot,'bin',computer('arch'),'libmwfftw3.3.dylib');
            % fftwheader = fullfile(matlabroot,'extern','include','fftw3.h');
            ipath = ['-I' fullfile(matlabroot,'extern','include')];

            % mex(ipath,'create_dct_plan_mex.cpp',fftwlibpath);
            mex(ipath,'fftw_plan_guru_r2r.cpp',fftwlibpath);
            mex(ipath,'execute_dct_plan_mex.cpp',fftwlibpath);
            mex(ipath,'execute_dct_plan_inout_mex.cpp',fftwlibpath);
        end

        function makeOMPMexFiles(fftwpath,omppath)
            arguments
                fftwpath = '/usr/local/opt/fftw'
                omppath = '/usr/local/opt/libomp'
            end
            lib_fftw = fullfile(fftwpath,'lib','libfftw3.dylib');
            lib_fftwomp = fullfile(fftwpath,'lib','libfftw3_omp.dylib');
            lib_omp = fullfile(omppath,'lib','libomp.dylib');
            inc_fftw = ['-I' fullfile(fftwpath,'include')];
            inc_omp = ['-I' fullfile(omppath,'include')];

            mex(inc_fftw,inc_omp,'create_dct_plan_mex.cpp',lib_fftw,lib_fftwomp,lib_omp);
            mex(inc_fftw,inc_omp,'execute_dct_plan_mex.cpp',lib_fftw,lib_fftwomp,lib_omp);
        end
    end
end