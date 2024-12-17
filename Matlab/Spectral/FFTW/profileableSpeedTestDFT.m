function profileableSpeedTestDFT
N = 1024;
NyNz = 256*256;
x = rand(N,NyNz);
nLoops = 10;

% DCT = WVTransformConstantStratification.CosineTransformForwardMatrix(N);
dft = RealToComplexTransform(size(x),dims=1,nCores=8,planner="measure");
% xout = zeros(size(x));

tic
for i=1:nLoops
    y1 = fft(x,N,1);
end
val = toc;
fprintf('%.2fs: matlab out-of-place, fft\n', val);

tic
for i=1:nLoops
    y2 = dft.transformForward(x);
end
val = toc;
fprintf('%.2fs: fftw out-of-place, r2c\n', val);

xout = complex(zeros(dft.complexSize));
tic
for i=1:nLoops
    xout = dft.transformForwardIntoArray(x,xout);
end
val = toc;
fprintf('%.2fs: fftw out-of-place, preallocated r2c\n', val);

tic
for i=1:nLoops
    x1back = ifft(y1,N,1,'symmetric');
end
val = toc;
fprintf('%.2fs: matlab out-of-place, ifft\n', val);

tic
for i=1:nLoops
    x2back = dft.transformBack(y2);
end
val = toc;
fprintf('%.2fs: fftw out-of-place, c2r\n', val);

xout2 = zeros(dft.realSize);
tic
for i=1:nLoops
    xout2 = dft.transformBackIntoArray(xout,xout2);
end
val = toc;
fprintf('%.2fs: fftw out-of-place, preallocated c2r\n', val);

N = 1024;
x = rand(N,N);
nLoops = 1000;
tic
for i=1:nLoops
    y1 = fft2(x);
end
val = toc;
fprintf('%.2fs: matlab out-of-place, fft2\n', val);

dft = RealToComplexTransform(size(x),dims=[1 2],nCores=8,planner="measure");

tic
for i=1:nLoops
    y2 = dft.transformForward(x);
end
val = toc;
fprintf('%.2fs: fftw out-of-place, r2c\n', val);

xout = complex(zeros(dft.complexSize));
tic
for i=1:nLoops
    xout = dft.transformForwardIntoArray(x,xout);
end
val = toc;
fprintf('%.2fs: fftw out-of-place, preallocated r2c\n', val);

end