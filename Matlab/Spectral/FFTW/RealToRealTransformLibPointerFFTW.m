classdef RealToRealTransformLibPointerFFTW < RealToRealTransform
    properties
        libname 
        libpath
        outp
    end

    methods
        function self = RealToRealTransformLibPointerFFTW(sz,options,newoptions)
            arguments
                sz
                options.dim double = 1
                options.transform char {mustBeMember(options.transform,["cosine","sine"])}
                options.planner char {mustBeMember(options.planner,["estimate","measure","patient","exhaustive"])} = "estimate"
                options.nCores = 1

                newoptions.libname = 'libmwfftw3'
                newoptions.libpath
            end

            scArgs = namedargs2cell(options);
            self@RealToRealTransform(sz,scArgs{:});

            self.libname = newoptions.libname;
            if isfield(newoptions,'libpath')
                self.libpath = newoptions.libpath;
            else
                self.libpath = fullfile(matlabroot,'bin',computer('arch'));
            end

            addpath(self.libpath);
            if not(libisloaded(self.libname))
                loadlibrary(self.libname,'fftw3.h')
            end

            in = zeros(sz);
            out = zeros(sz);
            calllib(self.libname,'fftw_plan_with_nthreads',options.nCores);
            self.plan = calllib(self.libname,'fftw_plan_guru_r2r',length(self.dims),self.dims,length(self.howmany_dims),self.howmany_dims,in,out,self.transformKind,uint32(self.planner));
            self.outp = libpointer('doublePtr',zeros(sz));
        end



        function f = transformBack(self,fbar)
            
            calllib(self.libname,'fftw_execute_r2r', self.plan, fbar, self.outp );
            f = self.outp.Value/2;
        end

        function fbar = transformForward(self,f)
            calllib(self.libname,'fftw_execute_r2r', self.plan, f, self.outp );
            fbar = self.scaleFactor*self.outp.Value;
        end
    end
end