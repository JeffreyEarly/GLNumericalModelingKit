function [f, xbar] = CosineTransformForward( tn, x )
% CosineTransformForward  Fast Discrete Cosine Transform (DST-I)
% 
% xbar is returned in the same units as x. This is the finite length
% definition of a Fourier transform.
%
% f is returned in units of cycles.
%
% The following relationship is satisfied:
%   [f, xbar] = CosineTransformForward(t,x);
%   S = T*(xbar .* conj(xbar));
%   x_sum = (1/T)*sum(x.*x)*dt;
%   S_sum = sum(S)*df;
%   x_sum == S_sum
%
% Using,
%   N=32; % total points
%   T=1.0; % total time length
%   t=T*(0:(N-1))'/N;

N = length(tn);

dctScratch = ifft(cat(1,x,x(N-1:-1:1)),2*N,1);
xbar = 2*real(dctScratch(2:N+1,:,:));

df = 1/(2*N*(tn(2)-tn(1)));
f = ((1:N)*df)';

end