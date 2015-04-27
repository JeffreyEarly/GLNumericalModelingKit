function uh = fct(u)
% Fast Chebyshev Transform 
%
%    uh : Discrete Chebyshev Transform Coefficients.
%    u  : Function values evaluted at Chebyshev Gauss Lobatto nodes
%         with nodes ordered increasingly x_i=[-1,...,1}
%         for i=1,2...,N
%
% By Allan P. Engsig-Karup, apek@imm.dtu.dk.

N  = length(u);
u  = ifft([u([1:N N-1:-1:2])]); % reverse ordering due to Matlab's fft
uh = ([u(1); 2*u(2:(N-1)); u(N)]);

if any(imag(uh))
	disp('Fast Chebyshev Transform returned imaginary values. Something went wrong!')
end

uh = real(uh);

return
