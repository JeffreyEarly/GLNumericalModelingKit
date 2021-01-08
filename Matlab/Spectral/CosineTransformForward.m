function [xbar, f] = CosineTransformForward( t, x, varargin )
% CosineTransformForward  Fast Discrete Cosine Transform (DCT-I)
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

if length(varargin) >= 1
    dim = varargin{1};
else
    dims = find(size(x) == length(t));
    if isempty(dims)
        error('Could not find a dimension with the same length as t.');
    elseif length(dims) == 1
        dim = dims;
    else
        error('You need to specifiy which dimension to differentiation, there are %d dimensions with the same length as t.',length(dims));
    end    
end

N = length(t);

% move the dim to the first dimension
% [x y dim] -> [dim x y]
x = shiftdim(x,dim-1);

switch ndims(x)
    case 2
        dctScratch = cat(1,x,x(N-1:-1:2,:));
        dctScratch = ifft(dctScratch,2*N-2,1);
        xbar = 2*(dctScratch(1:N,:));
        xbar(end,:) = xbar(end,:)/2;
    case 3
        dctScratch = cat(1,x,x(N-1:-1:2,:,:));
        dctScratch = ifft(dctScratch,2*N-2,1);
        xbar = 2*(dctScratch(1:N,:,:));
        xbar(end,:,:) = xbar(end,:,:)/2;
    otherwise
        error('Not yet implemented for more than 3 dimensions.')
end

% now move the dim back to where it was
% [dim x y]
xbar = shiftdim(xbar,ndims(xbar)-(dim-1));

if (max(abs(imag(xbar(:))))/max(abs(real(xbar(:)))) < 1e-14)
    xbar = real(xbar);
end

df = 1/(2*(N-1)*(t(2)-t(1)));
f = df*(0:(N-1))';

end