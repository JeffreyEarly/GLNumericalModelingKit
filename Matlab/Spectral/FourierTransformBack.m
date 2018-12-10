function [x, t] = FourierTransformBack( f, xbar, dim, varargin )

nT = length(xbar);
fourierFrequencyT = f(2)-f(1);
deltaT = 1/(nT*fourierFrequencyT);
% nyquistFrequencyT = 1/(2*deltaT);	% nyquist frequency

t = (0:nT-1)'*deltaT;

if length(varargin) == 1
x = ifft( xbar, [], dim,varargin{1} )*nT;
else
    x = ifft( xbar, [], dim )*nT;
end