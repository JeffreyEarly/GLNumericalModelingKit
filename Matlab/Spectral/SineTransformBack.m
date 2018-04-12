function [x, t] = SineTransformBack( f, xbar, includeEndpoints )
% SineTransformBack  Fast Discrete Inverse Sine Transform
% 

N = length(f);

dstScratch = 0.5*sqrt(-1)*cat(1,0,xbar,0,-xbar(N:-1:1));
x = fft(dstScratch,2*N+2,1);

x = real(x(1:N+2));
deltaT = 1/(2*(N+1)*(f(2)-f(1)));
t = (0:N+1)'*deltaT;



end