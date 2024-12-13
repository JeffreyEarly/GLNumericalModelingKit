classdef RealToRealTransformLibPointerFFTW < RealToRealTransform
    properties
        libname 
        libpath
        outp

        dims
        howmany_dims
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
                    self.dims(rank).n = alldims(iDim).n;
                    self.dims(rank).is = alldims(iDim).is;
                    self.dims(rank).os = alldims(iDim).os;
                    self.scaleFactor = 1/(self.dims(rank).n -1);
                else
                    howmany_rank = howmany_rank + 1;
                    self.howmany_dims(howmany_rank).n = alldims(iDim).n;
                    self.howmany_dims(howmany_rank).is = alldims(iDim).is;
                    self.howmany_dims(howmany_rank).os = alldims(iDim).os;
                end
            end

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