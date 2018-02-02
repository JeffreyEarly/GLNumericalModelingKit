%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Fourier transform
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N=32; % total points
T=4.0; % total time length
t=T*(0:(N-1))'/N;

% Create a simple function, x = 2 + 3*sin(4*pi*t);
% The lowest resolved frequency is one cycle per unit time
x = 2*ones(size(t)) + 3*sin(2*(2*pi)*t);

[f, xbar] = TransformForward(t,x,1);

S = T*(xbar .* conj(xbar));

dt = t(2)-t(1);
df = f(2)-f(1);

x_sum = (1/T)*sum(x.*x)*dt;
S_sum = sum(S)*df;

if abs(x_sum - S_sum) < 1e-7
   disp('Power matches')
else
   disp('Fail')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Sine transform
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create a simple function, x = sin(pi*t) + 3*sin(4*pi*t);
% The lowest resolved frequency is *half* a cycle per unit time
x =  sin(pi*t) + 3*sin(2*(2*pi)*t);

[f, xbar] = SineTransformForward(t,x);

S = T*(xbar .* conj(xbar));

df = f(2)-f(1);

x_sum = (1/T)*sum(x.*x)*dt;
S_sum = sum(S)*df;

if abs(x_sum - S_sum) < 1e-7
   disp('Power matches')
else
   disp('Fail')
end