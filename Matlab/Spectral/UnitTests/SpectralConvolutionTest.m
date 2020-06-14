%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Fourier transform
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Testing Fourier transform:\n')
N=16; % total points
T=4.0; % total time length
t=T*(0:(N-1))'/N;

% Create a simple function, x = 2 + 3*sin(4*pi*t);
% The lowest resolved frequency is one cycle per unit time
x = 1*ones(size(t)) + 3*sin(0.5*(2*pi)*t);

[xbar, f] = FourierTransformForward(t,x,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Convolution
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xbar = fftshift(xbar);
f = fftshift(f);

% f_n =       [-4 -3 -2 -1  0  1  2  3 ]
% flip(f_n) = [ 3  2  1  0 -1 -2 -3 -4 ]
% added together, we get an offset of -1 for all columns

% -8: [1] .* [1]
% -7: [1 2] .* [2 1]
% -6: [1 2 3] .* [3 2 1]
% ...
% -1: [1 2 3 4 5 6 7 8] .* [8 7 6 5 4 3 2 1]
%  0: [2 3 4 5 6 7 8] .* [8 7 6 5 4 3 2]
%  1: [3 4 5 6 7 8] .* [8 7 6 5 4 3]
% ...
%  7: []

% return;
% 
% xbar_f = flip(xbar);
xconv2N = zeros(2*N,1);
for i = 1:(2*N-1)
    indices = max(1,(i+1)-N):min(i,N);
    xconv2N(i) = sum ( xbar(indices) .* xbar(flip(indices)) );
end

xconvN = zeros(N,1);
for i = (N/2 + 1):(3*N/2)
    indices = max(1,(i+1)-N):min(i,N);
    xconvN(i-N/2) = sum ( xbar(indices) .* xbar(flip(indices)) );
end

xpseudo = fftshift(FourierTransformForward(t,x.*x,1));

% Now repeat, but without shifting, and only doing the first half
% [0  1  2  3 -4 -3 -2 -1]
% 0: [1 2 3 4] .* [1 8 7 6]
% 1: [1 2 3] .* [2 
% ooh, it's a mess with zero like this.
xbar = fftshift(xbar);
xconvNreal = zeros(N,1);


