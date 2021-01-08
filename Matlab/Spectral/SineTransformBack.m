function [x, t] = SineTransformBack( f, xbar, varargin )
% SineTransformBack  Fast Discrete Inverse Sine Transform
% 
% by default we're assuming that you're only doing resolved frequencies...
% so no nyquist and no zero.

if length(varargin) >= 1
    dim = varargin{1};
else
    dims = find(size(xbar) == length(f));
    if isempty(dims)
        error('Could not find a dimension with the same length as t.');
    elseif length(dims) == 1
        dim = dims;
    else
        error('You need to specifiy which dimension to differentiation, there are %d dimensions with the same length as t.',length(dims));
    end
end

N = length(f);

% move the dim to the first dimension
% [x y dim] -> [dim x y]
xbar = shiftdim(xbar,dim-1);

switch ndims(xbar)
    case 2
        dstScratch = 0.5*sqrt(-1)*cat(1,zeros(size(xbar(1,:))),xbar,zeros(size(xbar(1,:))),-xbar(N:-1:1,:));
        x = fft(dstScratch,2*N+2,1);
        x = (x(1:N+2,:)); % no longer taking the real part, to allow for complex functions
    case 3
        dstScratch = 0.5*sqrt(-1)*cat(1,zeros(size(xbar(1,:,:))),xbar,zeros(size(xbar(1,:,:))),-xbar(N:-1:1,:,:));
        x = fft(dstScratch,2*N+2,1);
        x = (x(1:N+2,:,:)); % no longer taking the real part, to allow for complex functions
    otherwise
        error('Not yet implemented for more than 3 dimensions.');
end

% now move the dim back to where it was
% [dim x y]
x = shiftdim(x,ndims(x)-(dim-1));

if (max(abs(imag(x(:))))/max(abs(real(x(:)))) < 1e-14)
    x = real(x);
end

deltaT = 1/(2*(N+1)*(f(2)-f(1)));
t = (0:N+1)'*deltaT;



end