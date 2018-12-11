function Dx = DiffFourier(t,x,varargin)
%DiffFourier Fourier transform derivative

if length(varargin) >= 1
    numDerivs = varargin{1};
else
    numDerivs = 1;
end

if length(varargin) >= 2
    dim = varargin{2};
else
    dims = find(size(x) == length(t));
    if isempty(dims)
        error('Could not find a dimension with the same length as t.');
    elseif length(dims) == 1
        dim = dims;
    else
        error('You need to specifiy which dimension to differentiation, there are (at least) two dimensions with the same length as t.');
    end
end

[xbar, f] = FourierTransformForward( t, x, dim );

% reshape the dimension so that we can multiply
newdims = ones(1,ndims(x));
newdims(dim) = length(f);
omega = reshape(2*pi*f,newdims);

Dx = FourierTransformBack(f, ((sqrt(-1)*omega).^numDerivs).*xbar, dim, 'symmetric');

end

