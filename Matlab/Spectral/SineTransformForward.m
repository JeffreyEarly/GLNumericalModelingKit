function [xbar,f] = SineTransformForward( t, x, varargin )
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

if length(varargin) >= 2
    endpointsIncluded = varargin{2};
else
    endpointsIncluded = 'both';
end

if length(varargin) >= 3
    shouldPerformSanityCheck = varargin{3};
else
    shouldPerformSanityCheck = 1;
end

% move the dim to the first dimension
% [x y dim] -> [dim x y]
x = shiftdim(x,dim-1);


if shouldPerformSanityCheck == 1
    reltol = 1e-13;
    abstol = 1e-13;
    mag = max(abs(x(:)));
    
    iszero = @(x) (x/mag) < reltol || x < abstol;
    
    switch ndims(x)
        case 2
            leftPoint = max(abs(x(1,:)));
            rightPoint = max(abs(x(end,:)));
        case 3
            leftPoint = max(abs(x(1,:)));
            rightPoint = max(abs(x(end,:)));
        otherwise
            error('Not yet implemented for more than 3 dimensions.')
    end

    if strcmp(endpointsIncluded,'both') == 1
        if ~iszero(leftPoint)|| ~iszero(rightPoint)
            fprintf('warning: by assumption both end points should be zero, but they are not: x(1)=%g, x(end)=%g.\n',x(1),x(end));
%             x(1) = 0; x(end) = 0;
        end
    elseif strcmp(endpointsIncluded,'left') == 1
        if ~iszero(leftPoint)
            fprintf('warning: by assumption the left point should be zero, but it is not: x(1)=%g.\n',x(1));
%             x(1) = 0;
        end
    elseif strcmp(endpointsIncluded,'right') == 1
        if ~iszero(rightPoint)
            fprintf('warning: by assumption the right point should be zero, but it is not: x(end)=%g.\n',x(end));
%             x(end) = 0;
        end
    end
end


if strcmp(endpointsIncluded,'both') == 1
    N = length(t)-1;
    switch ndims(x)
        case 2
            dstScratch = cat(1,x,-x(N:-1:2,:)); % 0, a, b, c, 0, -c, -b, -a
        case 3
            dstScratch = cat(1,x,-x(N:-1:2,:,:)); % 0, a, b, c, 0, -c, -b, -a
        otherwise
            error('Not yet implemented for more than 3 dimensions.')
    end
elseif strcmp(endpointsIncluded,'left') == 1
    N = length(t);
    switch ndims(x)
        case 2
            dstScratch = cat(1,x,0,-x(N:-1:2,:));  % 0, a, b, c, 0, -c, -b, -a
        case 3
            dstScratch = cat(1,x,0,-x(N:-1:2,:,:));  % 0, a, b, c, 0, -c, -b, -a
        otherwise
            error('Not yet implemented for more than 3 dimensions.')
    end
elseif strcmp(endpointsIncluded,'right') == 1    
    N = length(t);
    switch ndims(x)
        case 2
            dstScratch = cat(1,0,x,-x(N-1:-1:1,:));  % 0, a, b, c, 0, -c, -b, -a
        case 3
            dstScratch = cat(1,0,x,-x(N-1:-1:1,:,:));  % 0, a, b, c, 0, -c, -b, -a
        otherwise
            error('Not yet implemented for more than 3 dimensions.')
    end
elseif strcmp(endpointsIncluded,'none') == 1
    N = length(t)+1;
    switch ndims(x)
        case 2
            dstScratch = cat(1,0,x,0,-x(N-1:-1:1,:));   % 0, a, b, c, 0, -c, -b, -a
        case 3
            dstScratch = cat(1,0,x,0,-x(N-1:-1:1,:,:));   % 0, a, b, c, 0, -c, -b, -a
        otherwise
            error('Not yet implemented for more than 3 dimensions.')
    end
end

dstScratch = ifft(dstScratch,2*N,1);

switch ndims(x)
    case 2
        xbar = -2*sqrt(-1)*(dstScratch(2:N,:));
    case 3
        xbar = -2*sqrt(-1)*(dstScratch(2:N,:,:));
    otherwise
        error('Not yet implemented for more than 3 dimensions.')
end
% now move the dim back to where it was
% [dim x y]
xbar = shiftdim(xbar,ndims(xbar)-(dim-1));

if (max(abs(imag(xbar(:))))/max(abs(real(xbar(:)))) < 1e-14)
    xbar = real(xbar);
end

df = 1/(2*N*(t(2)-t(1)));
f = ((1:N-1)*df)';

end