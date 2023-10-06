%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Fourier transform
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Testing Fourier transform:\n')
N=32; % total points
T=1.0; % total time length
t=T*(0:(N-1))'/N;

% Create a simple function, x = 2 + 3*sin(4*pi*t);
% The lowest resolved frequency is one cycle per unit time
x = 1*ones(size(t)) + 3*cos(16*(2*pi)*t);
y = 5*ones(size(t)) + 2*cos(16*(2*pi)*t);

xbar = FourierTransformForward(t,x);
ybar = FourierTransformForward(t,y);

c=x+sqrt(-1)*y;
[cbar, f] = FourierTransformForward(t,c);

xbar2 = zeros(N/2+1,1);
ybar2 = zeros(N/2+1,1);
xbar2(1) = 2*real(cbar(1));
ybar2(1) = 2*imag(cbar(1));
for i=1:N/2
   xbar2(i+1) = cbar(i+1)+conj(cbar(N-i+1));
   ybar2(i+1) = sqrt(-1)*(-cbar(i+1)+conj(cbar(N-i+1)));
end


return

S = T*(cbar .* conj(cbar));

dt = t(2)-t(1);
df = f(2)-f(1);

x_sum = (1/T)*sum(c.*conj(c))*dt;
S_sum = sum(S)*df;

if abs(x_sum - S_sum) < 1e-7
   fprintf('\tPower matches.\n')
else
   fprintf('\tPower FAIL\n')
end

[c_back, t_back] = FourierTransformBack(f,cbar,1);
if (max(abs(c-c_back)) < 1e-7) && (max(abs(t-t_back)) < 1e-7)
   fprintf('\tTransformed backed.\n')
else
   fprintf('\tInverse transform FAIL.\n')
end


