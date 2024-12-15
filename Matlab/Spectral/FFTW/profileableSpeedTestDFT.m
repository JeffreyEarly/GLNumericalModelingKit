function profileableSpeedTestDFT
N = 512;
NyNz = 256*256;
x = rand(N,NyNz);
nLoops = 10;

% DCT = WVTransformConstantStratification.CosineTransformForwardMatrix(N);
dft = RealToComplexTransform(size(x),dims=1,nCores=8); % ,planner="measure"
% xout = zeros(size(x));

tic
for i=1:nLoops
    y1 = fft(x,N,1);
end
val = toc;
fprintf('%.2f: matlab out-of-place, fft\n', val);

tic
for i=1:nLoops
    y2 = dft.transformForward(x);
end
val = toc;
fprintf('%.2f: fftw out-of-place, r2c\n', val);

xout = complex(zeros(dft.complexSize));
tic
for i=1:nLoops
    xout = dft.transformForwardIntoArray(x,xout);
end
val = toc;
fprintf('%.2f: fftw out-of-place, preallocated r2c\n', val);

% tic
% for i=1:nLoops
%     xout = dct.transformForwardIntoArray(x,xout);
% end
% val = toc;
% fprintf('%.2f: out-of-place preallocated, fast DFT\n', val);


% tic
% for i=1:nLoops
%     % x = dct.transformForward(x);
%     % x = execute_dct_plan_mex(dct.plan,x);
%     % xout = execute_dct_plan_inout_mex(dct.plan,x,xout);
%     % xout = dct.transformForwardIntoArray(x,xout);
%     % x = dct.transformForward(x);
%     % x = applyMatrixDCT(x,DCT);
%     % y = unaryOperation(x);
%     % y = removeNegativeNumbers(x);
% end
% toc

end