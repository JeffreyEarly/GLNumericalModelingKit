function [ S_blur ] = BlurSpectrum( omega, S )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
N=length(omega);
MF = floor(N/2)+1;
SpecPermute = ifftshift(S); %[S(MF:N); S(1:MF-1)]; % Permute Spectrum for ifft
Autocov = ifft(SpecPermute); % Autocovariance of Spectrum
% TriangleKernel = ifftshift((1-([ceil(N/2):-1:1 0:MF-2]')/MF));

TriangleKernel = (1-([0:MF-1 ceil(N/2)-1:-1:1]/MF))'; % Triange kernel (works for odd and even N)
S_blur = fftshift(real(fft(TriangleKernel.*Autocov))); % Blurred Spectrum
end

