function [xbar,f] = SineTransformForward( t, x, endpointsIncluded )
% SineTransformForward  Fast Discrete Sine Transform (DST-I)
% 
% xbar is returned in the same units as x. This is the finite length
% definition of a Fourier transform.
%
% f is returned in units of cycles, and only contains the potentially
% nonzero frequency i.e., the zero and nyquist are not included, because
% they must be zero.
%
% The following relationship is satisfied:
%   [f, xbar] = SineTransformForward(t,x);
%   S = T*(xbar .* conj(xbar));
%   x_sum = (1/T)*sum(x.*x)*dt;
%   S_sum = sum(S)*df;
%   x_sum == S_sum
%
% Using,
%   N=33; % total points
%   T=1.0; % total time length
%   t=T*(0:(N-1))'/(N-1);
% This definition of t *includes* the end points which *must* be zero.
%
%


if nargin < 3
    endpointsIncluded = 'both'; 
end

eps = 1e-14;
mag = max(abs(x));
if strcmp(endpointsIncluded,'both') == 1
    if abs(x(1))/mag > eps || abs(x(end))/mag > eps
        fprintf('warning: by assumption both end points should be zero, but they are not: x(1)=%g, x(end)=%g.\n',x(1),x(end));
        x(1) = 0; x(end) = 0;
    end
    
    N = length(t)-1;
    dstScratch = cat(1,x,-x(N:-1:2)); % 0, a, b, c, 0, -c, -b, -a
elseif strcmp(endpointsIncluded,'left') == 1
    if abs(x(1))/mag > eps
        fprintf('warning: by assumption the left point should be zero, but it is not: x(1)=%g.\n',x(1));
        x(1) = 0;
    end
    
    N = length(t);
    dstScratch = cat(1,x,0,-x(N:-1:2));  % 0, a, b, c, 0, -c, -b, -a
elseif strcmp(endpointsIncluded,'right') == 1
    if abs(x(end))/mag > eps
        fprintf('warning: by assumption the right point should be zero, but it is not: x(end)=%g.\n',x(end));
        x(end) = 0;
    end
    
    N = length(t);
    dstScratch = cat(1,0,x,-x(N-1:-1:1));  % 0, a, b, c, 0, -c, -b, -a
elseif strcmp(endpointsIncluded,'none') == 1
    N = length(t)+1;
    dstScratch = cat(1,0,x,0,-x(N-1:-1:1));   % 0, a, b, c, 0, -c, -b, -a
end

dstScratch = ifft(dstScratch,2*N,1);
xbar = 2*imag(dstScratch(2:N,:,:));

df = 1/(2*N*(t(2)-t(1)));
f = ((1:N-1)*df)';

end