function u = ifct(uh)
% Inverse Fast Chebyshev Transform 
%
%    uh : Discrete Chebyshev Transform Coefficients.
%    u  : Function values evaluted at Chebyshev Gauss Lobatto nodes
%         with nodes ordered increasingly x_i=[-1,...,1}
%         for i=1,2...,N
%
% By Allan P. Engsig-Karup, apek@imm.dtu.dk.

N = length(uh);
u=fft([uh(1); [uh(2:N-1); uh(N)*2; uh(N-1:-1:2)]*0.5]); 
u=(u(1:N)); % reverse ordering due to Matlab's fft.

if any(imag(u)./abs(u) > 1e-6)
	disp('Inverse Fast Chebyshev Transform returned imaginary values. Something went wrong!')
end

u = real(u);

return
