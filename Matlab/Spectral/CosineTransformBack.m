function [x, t] = CosineTransformBack( f, xbar )
% CosineTransformBack  Fast Discrete Inverse Cosine Transform
% 

N = length(f);

dctScratch = cat(1,0.5*xbar(1:N-1),xbar(N),0.5*xbar(N-1:-1:2));
x = fft(dctScratch,2*N-2,1);

x = real(x(1:N));
deltaT = 1/(2*(N-1)*(f(2)-f(1)));
t = (0:N-1)'*deltaT;
end