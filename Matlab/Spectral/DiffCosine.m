function Dx = DiffCosine(t,x,varargin)
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

[xbar, f] = CosineTransformForward( t, x, dim );

% reshape the dimension so that we can multiply
newdims = ones(1,ndims(x));
newdims(dim) = length(f);
f = reshape(f,newdims);

Dxbar = ((sqrt(-1)*2*pi*f).^numDerivs).*xbar;

if mod(numDerivs,2) == 0
    Dx = CosineTransformBack(f, Dxbar, dim);
else
    Dx = SineTransformBack(f, circshift(Dxbar,-1,dim), dim);
end



end

