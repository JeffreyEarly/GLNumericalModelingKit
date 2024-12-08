% mex '-LC:\fftw-3.3.4-dll32' -llibfftw3-3.lib test.c

fftwlib = fullfile(matlabroot,'bin',computer('arch'),'libmwfftw3.3.dylib');
fftwheader = fullfile(matlabroot,'extern','include','fftw3.h');
ipath = ['-I' fullfile(matlabroot,'extern','include')];

% mex(ipath,'dctFFTW.c',fftwlib)
%%
mex(ipath,'create_dct_plan_mex.cpp',fftwlib)

%%
mex(ipath,'execute_dct_plan_mex.cpp',fftwlib)

%%
N = 257;
NyNz = 10*256*256;
x = rand(N,NyNz);
nLoops = 20;
sz = size(x);
options.dim = 1;

DCT = WVTransformConstantStratification.CosineTransformForwardMatrix(N);

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

plan = create_dct_plan_mex(rank, dims, howmany_rank, howmany_dims);

%%
y = zeros(size(x));
tic
for i=1:nLoops
    execute_dct_plan_mex(plan,x);
end
toc

%%
y = zeros(size(x));
tic
for i=1:nLoops
    y = DCT*x;
end
toc