function [x, t] = CosineTransformBack( f, xbar, varargin )
% CosineTransformBack  Fast Discrete Inverse Cosine Transform
% 

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
        dctScratch = cat(1,0.5*xbar(1:N-1,:),xbar(N,:),0.5*xbar(N-1:-1:2,:));
%         x = fft(dctScratch,2*N-2,1);
        x = ifft(dctScratch,2*N-2,1)*(2*N-2);
        x = (x(1:N,:)); % no longer taking the real part, to allow for complex functions
    case 3
        dctScratch = cat(1,0.5*xbar(1:N-1,:,:),xbar(N,:,:),0.5*xbar(N-1:-1:2,:,:));
%         x = fft(dctScratch,2*N-2,1);
        x = ifft(dctScratch,2*N-2,1)*(2*N-2);
        x = (x(1:N,:,:)); % no longer taking the real part, to allow for complex functions
    otherwise
        error('Not yet implemented for more than 3 dimensions.');
end

% now move the dim back to where it was
% [dim x y]
x = shiftdim(x,ndims(x)-(dim-1));

if (max(abs(imag(x(:))))/max(abs(real(x(:)))) < 1e-14)
    x = real(x);
end

deltaT = 1/(2*(N-1)*(f(2)-f(1)));
t = (0:N-1)'*deltaT;
end