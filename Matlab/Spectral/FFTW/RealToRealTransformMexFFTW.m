classdef RealToRealTransformMexFFTW < RealToRealTransform
    methods
        function self = RealToRealTransformMexFFTW(sz,options)
            arguments
                sz
                options.dim double = 1
                options.transform char {mustBeMember(options.transform,["cosine","sine"])}
                options.planner char {mustBeMember(options.planner,["estimate","measure","patient","exhaustive"])} = "estimate"
                options.nCores = 1
            end

            scArgs = namedargs2cell(options);
            self@RealToRealTransform(sz,scArgs{:});

            self.plan = create_dct_plan_mex(self.dims, self.howmany_dims, options.nCores, self.transformKind, self.planner);
        end

        function f = transformBack(self,fbar)     
            f = execute_dct_plan_mex(self.plan,fbar)/2;
        end

        % function fbar = transformForward(self,f)
        %     fbar = self.scaleFactor*execute_dct_plan_mex(self.plan,f);
        % end
        function f = transformForward(self,f)
            % fout = zeros(size(f));
            % fout = self.scaleFactor*execute_dct_plan_mex(self.plan,f,fout);
            f = self.scaleFactor*execute_dct_plan_mex(self.plan,f);
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

            mex(ipath,'create_dct_plan_mex.cpp',fftwlibpath);
            mex(ipath,'execute_dct_plan_mex.cpp',fftwlibpath);
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