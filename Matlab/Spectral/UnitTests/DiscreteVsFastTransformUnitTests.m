%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Cosine transform (Real)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = 33;
T=4.0; % total time length
t=T*(0:(N-1))'/(N-1);

x = randn(size(t));

[xbar, f] = CosineTransformForward(t,x);
DCT = CosineTransformForwardMatrix(N);
xbar2 = DCT*x;

epsilon = 1 - xbar./xbar2;
if (max(abs(epsilon)) < 1e-7)
   fprintf('\tForward DCT matches.\n')
else
   fprintf('\tForward DCT FAIL.\n')
end

[x_back, t_back] = CosineTransformBack(f,xbar);
iDCT = CosineTransformBackMatrix(N);
x_back2 = iDCT*xbar2;

epsilon = 1 - x_back./x_back2;
if (max(abs(epsilon)) < 1e-7)
   fprintf('\tInverse DCT matches.\n')
else
   fprintf('\tInverseDCT FAIL.\n')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Sine transform (Real)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = 33;
T=4.0; % total time length
t=T*(0:(N-1))'/(N-1);

x = randn(size(t));
x(1) = 0; x(end) = 0;

[xbar, f] = SineTransformForward(t,x);
DST = SineTransformForwardMatrix(N);
xbar2 = DST*x;

epsilon = 1 - xbar./xbar2;
if (max(abs(epsilon)) < 1e-7)
   fprintf('\tForward DST matches.\n')
else
   fprintf('\tForward DST FAIL.\n')
end

[x_back, t_back] = SineTransformBack(f,xbar);
iDST = SineTransformBackMatrix(N);
x_back2 = iDST*xbar2;

epsilon = x_back - x_back2;
if (max(abs(epsilon)) < 1e-7)
   fprintf('\tInverse DST matches.\n')
else
   fprintf('\tInverse DST FAIL.\n')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Fourier transform (Real)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N=32; % total points
T=4.0; % total time length
t=T*(0:(N-1))'/N;

x = randn(size(t));

[xbar, f] = FourierTransformForward(t,x);
DFT = FourierTransformForwardMatrix(N);
rDFT = FourierTransformForwardMatrixReal(N);
xbar2 = DFT*x;
xbar3 = rDFT*x;

epsilon = xbar - xbar2;
if (max(abs(epsilon)) < 1e-7)
   fprintf('\tForward DFT matches.\n')
else
   fprintf('\tForward DFT FAIL.\n')
end
epsilon = xbar(1:(N/2+1)) - xbar3;
if (max(abs(epsilon)) < 1e-7)
   fprintf('\tForward real DFT matches.\n')
else
   fprintf('\tForward real DFT FAIL.\n')
end

[x_back, t_back] = FourierTransformBack(f,xbar,1);
iDFT = FourierTransformBackMatrix(N);
irDFT = FourierTransformBackMatrixReal(N);
x_back2 = iDFT*xbar2;

x_back3 = real(irDFT*xbar3);

epsilon = x_back - x_back2;
if (max(abs(epsilon)) < 1e-7)
   fprintf('\tInverse DFT matches.\n')
else
   fprintf('\tInverse DFT FAIL.\n')
end

epsilon = x_back - x_back3;
if (max(abs(epsilon)) < 1e-7)
   fprintf('\tInverse real DFT matches.\n')
else
   fprintf('\tInverse real DFT FAIL.\n')
end