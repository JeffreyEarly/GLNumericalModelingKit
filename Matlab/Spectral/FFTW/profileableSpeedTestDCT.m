function profileableSpeedTestDCT
N = 513;
NyNz = 256*256;
x = rand(N,NyNz);
nLoops = 10;

DCT = WVTransformConstantStratification.CosineTransformForwardMatrix(N);
dct = RealToRealTransformMexFFTW(size(x),dim=1,transform="cosine",nCores=8);
xout = zeros(size(x));

tic
for i=1:nLoops
    x = DCT*x;
end
val = toc;
fprintf('%.2fs: in-place, matrix DCT\n', val);

tic
for i=1:nLoops
    y = dct.transformForward(x);
end
val = toc;
fprintf('%.2fs: out-of-place, fast DCT\n', val);

tic
for i=1:nLoops
    x = dct.transformForward(x);
end
val = toc;
fprintf('%.2fs: in-place, fast DCT\n', val);

tic
for i=1:nLoops
    xout = dct.transformForwardIntoArray(x,xout);
end
val = toc;
fprintf('%.2fs: out-of-place preallocated, fast DCT\n', val);

end