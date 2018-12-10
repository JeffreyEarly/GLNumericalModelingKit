function Dx = DiffSine(t,x,varargin)
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

% move the dim to the first dimension
% [x y dim] -> [dim x y]
x = shiftdim(x,dim-1);

[xbar, f] = SineTransformForward( t, x, 1 );

% reshape the dimension so that we can multiply
f = reshape(f,[],1);

% sine goes to, [1,-1,-1,1] for numDerivs = [1,2,3,4]
thesign = [1,-1,-1,1];
Dxbar = thesign(mod(numDerivs-1,4)+1)*((2*pi*f).^numDerivs).*xbar;

if mod(numDerivs,2) == 0
    Dx = SineTransformBack(f, Dxbar, 1);
else
    df = f(2)-f(1);
    f = cat(1,0,f,f(end)+df);
    
    
    switch ndims(Dxbar)
        case 2
            Dxbar = cat(1,zeros(size(Dxbar(1,:))),Dxbar,zeros(size(Dxbar(1,:))));
        case 3
            Dxbar = cat(1,zeros(size(Dxbar(1,:,:))),Dxbar,zeros(size(Dxbar(1,:,:))));
        otherwise
            error('Not yet implemented for more than 3 dimensions.');
    end
    
    Dx = CosineTransformBack(f, Dxbar, 1);
end

% now move the dim back to where it was
% [dim x y]
Dx = shiftdim(Dx,ndims(Dx)-(dim-1));

end

