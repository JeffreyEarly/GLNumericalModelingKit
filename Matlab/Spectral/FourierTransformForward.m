% t = 0:1/64:(1-1/64);
% a = sin(2*pi*3*t)';
% [fn, Yn] = TransformForward( t, a, 2);
% figure, plot( fn, real(Yn), 'blue'), hold on, plot(fn, imag(Yn), 'red')
% [t2, b] = TransformBack( fn, Yn, 1);
% figure, plot(t, a, 'blue'), hold on, plot(t,real(b)+0.1,'red')
% c = 2*pi*sqrt(-1)*fn.*Yn;
% [t2, b] = TransformBack( fn, c, 2);
% figure, plot(t, a, 'blue'), hold on, plot(t,real(b),'red')

function [xbar, f] = FourierTransformForward( t, x, dim )
% FourierTransformForward Fast Fourier Transform
%
% xbar is returned in the same units as x. This is the finite length
% definition of a Fourier transform.
%
% f is returned in units of cycles.
%
% Using,
%   N=32; % total points
%   T=1.0; % total time length
%   t=T*(0:(N-1))'/N;
% we have defined t such that the signal is assumed periodic, so that the
% following relationship is satisfied:
%   [f, xbar] = FourierTransformForward(t,x);
%   S = T*(xbar .* conj(xbar));
%   x_sum = (1/T)*sum(x.*x)*dt;
%   S_sum = sum(S)*df;
%   x_sum == S_sum
%

nT = length(t);
f = FFTFrequenciesFromTimeSeries( t );
xbar = fft( x, nT, dim )/nT;
