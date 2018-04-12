function [xbar, f] = CosineTransformForward( t, x )
% CosineTransformForward  Fast Discrete Cosine Transform (DST-I)
% 
% xbar is returned in the same units as x. This is the finite length
% definition of a Fourier transform.
%
% f is returned in units of cycles.
%
% The cosine series would have the following sum,
%   x(t) = xbar(1)/2 + sum( xbar(i)*cos(i*pi/T) )
% So note that the first coefficient double the average of the function.
%
% From this, Parseval's theorem is,
%   x_sum = (1/T)*(sum(x(2:end-1).*x(2:end-1))*dt+x(1)*x(1)*dt/2 + x(end)*x(end)*dt/2);
%   S_sum = ( S(1)/2 + sum(S(2:end-1)) + 2*S(end))*df;
%
%   Where we've taken care to integrate only to the endpoints (and not
%   beyond) in the x_sum. The S_sum includes a correction for the Nyquist.

N = length(t);

dctScratch = cat(1,x,x(N-1:-1:2));
dctScratch = ifft(dctScratch,2*N-2,1);
xbar = 2*real(dctScratch(1:N));

xbar(end) = xbar(end)/2;

df = 1/(2*(N-1)*(t(2)-t(1)));
f = df*(0:(N-1))';

end