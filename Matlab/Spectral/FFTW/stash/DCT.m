if not(libisloaded('libmwfftw3'))
    addpath(fullfile(matlabroot,'bin','maca64'))
    loadlibrary('libmwfftw3','fftw3.h')
end

calllib('libmwfftw3','fftw_plan_with_nthreads',2)

iodim.n = 17;
iodim.is = 1;
iodim.os = 1;

loopdim.n = 0;
loopdim.is = 1;
loopdim.os = 1;

rank = 1;
transformKind = 3;
planner = bitshift(1,6);
in = zeros(iodim.n,1);
out = zeros(iodim.n,1);
inp = libpointer('doublePtr',in);
outp = libpointer('doublePtr',out);

% fftw_plan fftw_plan_guru_r2r(int rank, const fftw_iodim *dims,
%                              int howmany_rank,
%                              const fftw_iodim *howmany_dims,
%                              double *in, double *out,
%                              const fftw_r2r_kind *kind,
%                              unsigned flags);

plan = calllib('libmwfftw3','fftw_plan_guru_r2r',rank,iodim,0,loopdim,in,out,transformKind,planner);

%%
in = zeros(iodim.n,1);
in(2) = 1;
inp = libpointer('doublePtr',in);
calllib('libmwfftw3','fftw_execute_r2r', plan, inp, outp );
outp.Value/2

D = WVTransformConstantStratification.CosineTransformBackMatrix(iodim.n);
D*in