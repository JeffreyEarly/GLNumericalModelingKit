function profileableSpeedTest
N = 513;
NyNz = 256*256;
x = rand(N,NyNz);
nLoops = 10;

% DCT = WVTransformConstantStratification.CosineTransformForwardMatrix(N);
dct = RealToRealTransformMexFFTW(size(x),dim=1,transform="cosine",nCores=8);
xout = zeros(size(x));

tic
for i=1:nLoops
    % x = dct.transformForward(x);
    % x = execute_dct_plan_mex(dct.plan,x);
    % xout = execute_dct_plan_inout_mex(dct.plan,x,xout);
    % xout = dct.transformForwardIntoArray(x,xout);
    x = dct.transformForward(x);
    % x = applyMatrixDCT(x,DCT);
    % y = unaryOperation(x);
    % y = removeNegativeNumbers(x);
end
toc

end