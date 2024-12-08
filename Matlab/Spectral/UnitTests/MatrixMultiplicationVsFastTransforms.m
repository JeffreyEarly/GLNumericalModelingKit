%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Cosine transform
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = 257;
t=(0:(N-1))'/(N-1);

DCT = WVTransformConstantStratification.CosineTransformForwardMatrix(N);
iDCT = WVTransformConstantStratification.CosineTransformBackMatrix(N);

NyNz = 256*256;
x = rand(N,NyNz);
nLoops = 20;
% dct = RealToRealTransformLibPointerFFTW(size(x),dim=1,transform="cosine",planner="exhaustive",nCores=8);
% dct = RealToRealTransformLibPointerFFTW(size(x),dim=1,transform="cosine",planner="measure",nCores=8,libname='libfftw3_omp',libpath='/opt/homebrew/opt/fftw/lib');

% RealToRealTransformMexFFTW.makeMexFiles();
% RealToRealTransformMexFFTW.makeOMPMexFiles();
% dct = RealToRealTransformMexFFTW(size(x),dim=1,transform="cosine",planner="exhaustive",nCores=8);

%%
fprintf('%d point DCT/iDCT with matrix transform:\n\t', N)
tic
for i=1:nLoops
    xbar = DCT*x;
    xback = iDCT*xbar;
end
toc

% fprintf('%d point DCT/iDCT via FFT:\n\t', N)
% tic
% for i=1:nLoops
%     xbar2 = CosineTransformForward(t,x);
%     xback2 = CosineTransformBack(t,xbar2);
% end
% toc

fprintf('%d point DCT/iDCT via FFTW:\n\t', N)

%%
% x_ptr = libpointer('doublePtr',x);
% x_bar2_ptr = libpointer('doublePtr',zeros(size(x)));
% xback2_ptr = libpointer('doublePtr',zeros(size(x)));
% a = zeros(size(x));
% b = zeros(size(x));

N = 257;
NyNz = 256*256;
x = rand(N,NyNz);
nLoops = 20;

y = zeros(size(x));
tic
% profile on
for i=1:nLoops
    % calllib('libfftw3_omp','fftw_execute_r2r', dct.plan, x_ptr, x_bar2_ptr );
    % a = dct.scaleFactor*x_bar2_ptr.Value; % this .Value kills performance. It's copying data.
    % calllib('libfftw3_omp','fftw_execute_r2r', dct.plan, x_bar2_ptr, xback2_ptr );
    % b = xback2_ptr.Value/2;
    % x = dct.transformForward(x);
    % x = dct.transformBack(x);
    % x(x<0) = 0;
    x = unaryOperation(x);
    % y = simple_add(y,x);
    % x = execute_dct_plan_mex(dct.plan,x);
    % x = execute_dct_plan_mex(dct.plan,x);
end
% profile viewer
toc

return

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% %   Fourier transform
% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% N = 512;
% t=(0:(N-1))'/(N);
% 
% DFT = FourierTransformForwardMatrix(N);
% iDFT = FourierTransformBackMatrix(N);
% 
% rDFT = FourierTransformForwardMatrixReal(N);
% irDFT = FourierTransformBackMatrixReal(N);
% irDFTrp = real(irDFT);
% irDFTip = -imag(irDFT);
% 
% NyNz = 128*128;
% x = rand(N,NyNz);
% nLoops = 10;
% 
% fprintf('\n%d point DFT with matrix transform:\n\t', N)
% tic
% for i=1:nLoops
%      xbar = DFT*x;
% end
% toc
% 
% fprintf('%d point DFT via FFT:\n\t', N)
% tic
% for i=1:nLoops
%     xbar2 = FourierTransformForward(t,x);
% end
% toc
% 
% fprintf('%d point *real* DFT with matrix transform:\n\t', N)
% tic
% for i=1:nLoops
%     xbar3 = rDFT*x;
% end
% toc
% 
% 
% 
% fprintf('\n%d point iDFT with matrix transform:\n\t', N)
% xbar = DFT*x;
% tic
% for i=1:nLoops
% %     xbar = DFT*x;
%     xback = iDFT*xbar;
% end
% toc
% 
% fprintf('%d point iDFT via FFT:\n\t', N)
% xbar2 = FourierTransformForward(t,x);
% tic
% for i=1:nLoops
% %     xbar2 = FourierTransformForward(t,x);
%     xback2 = ifft(xbar2,'symmetric')*N;
% end
% toc
% 
% fprintf('%d point *real* iDFT with matrix transform:\n\t', N)
% xbar3 = rDFT*x;
% tic
% for i=1:nLoops
% %     xbar3 = rDFT*x;
%     xback3 = irDFTrp*real(xbar3) + irDFTip*imag(xbar3);
% end
% toc
% 
% fprintf('%d point iDFT with forced Hermitian:\n\t', N)
% xbar3 = rDFT*x;
% tic
% for i=1:nLoops
%     xback3 = ifft(xbar3,N,'symmetric')*N;
% end
% toc
% 

shouldAntiAlias = 1;

Ns = 2.^(6:10);
method = zeros(3,length(Ns));

for iN=1:length(Ns)
    N = Ns(iN);
    t=(0:(N-1))'/(N);

    DFT = FourierTransformForwardMatrix(N);
    iDFT = FourierTransformBackMatrix(N);
%     if (shouldAntiAlias == 1)
%         DFT = DFT(1:(N/2+1),:);
%         iDFT = iDFT(:,1:(N/2+1));
%     end

    rDFT = FourierTransformForwardMatrixReal(N);
    irDFT = FourierTransformBackMatrixReal(N);
    if (shouldAntiAlias == 1)
        K = floor(2*(N/2+1)/3);
        rDFT = rDFT(1:K,:);
        irDFT = irDFT(:,1:K);
    end
    irDFTrp = real(irDFT);
    irDFTip = -imag(irDFT);

    NyNz = 128*128;
    x = rand(N,NyNz);
    nLoops = 10;



    fprintf('\n\n%d point DFT/iDFT via FFT:\n\t', N)
    tic
    for i=1:nLoops
        xbar2 = FourierTransformForward(t,x);
        xback2 = ifft(xbar2,'symmetric')*N;
    end
    method(1,iN)=toc;

    fprintf('%d point *real* DFT/iDFT with matrix transform:\n\t', N)
    tic
    for i=1:nLoops
        xbar3 = rDFT*x;
        xback3 = irDFTrp*real(xbar3) + irDFTip*imag(xbar3);
    end
    method(2,iN)=toc;

    fprintf('%d point real DFT, iDFT with forced Hermitian:\n\t', N)
    tic
    for i=1:nLoops
        xbar3 = rDFT*x;
        xback3 = ifft(xbar3,N,'symmetric')*N;
    end
    method(3,iN)=toc;
end
figure, plot(log2(Ns),method.')
legend('FFT','MM','Hybrid')
