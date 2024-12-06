% mex '-LC:\fftw-3.3.4-dll32' -llibfftw3-3.lib test.c

fftwlib = fullfile(matlabroot,'bin',computer('arch'),'libmwfftw3.3.dylib');
fftwheader = fullfile(matlabroot,'extern','include','fftw3.h');
ipath = ['-I' fullfile(matlabroot,'extern','include')];

mex(ipath,'dctFFTW.c',fftwlib)