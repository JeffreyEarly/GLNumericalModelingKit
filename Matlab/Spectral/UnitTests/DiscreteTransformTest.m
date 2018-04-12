%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Fourier transform
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Testing Fourier transform:\n')
N=32; % total points
T=4.0; % total time length
t=T*(0:(N-1))'/N;

% Create a simple function, x = 2 + 3*sin(4*pi*t);
% The lowest resolved frequency is one cycle per unit time
x = 1*ones(size(t)) + 3*sin(2*(2*pi)*t);

[xbar, f] = FourierTransformForward(t,x,1);

S = T*(xbar .* conj(xbar));

dt = t(2)-t(1);
df = f(2)-f(1);

x_sum = (1/T)*sum(x.*x)*dt;
S_sum = sum(S)*df;

if abs(x_sum - S_sum) < 1e-7
   fprintf('\tPower matches.\n')
else
   fprintf('\tPower fail\n')
end

[x_back, t_back] = FourierTransformBack(f,xbar,1);
if (max(abs(x-x_back)) < 1e-7) && (max(abs(t-t_back)) < 1e-7)
   fprintf('\tTransformed backed.\n')
else
   fprintf('\tInverse transform fail.\n')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Sine transform
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Testing sine transform:\n')

N=33; % total points
T=4.0; % total time length
t=T*(0:(N-1))'/(N-1);

% Create a simple function, x = sin(pi*t) + 3*sin(4*pi*t);
% The lowest resolved frequency is *half* a cycle per unit time
x =  sin(pi*t/T) + 2*sin(8*(2*pi)*t/T) + 3*sin(15.5*(2*pi)*t/T);

[xbar,f] = SineTransformForward(t,x);

S = T*(xbar .* conj(xbar));

df = f(2)-f(1);

x_sum = (1/T)*sum(x.*x)*dt;
S_sum = sum(S)*df;

if abs(x_sum - S_sum) < 1e-7
   fprintf('\tPower matches\n')
else
   fprintf('\tPower fail\n')
end

[x_back, t_back] = SineTransformBack(f,xbar);
if (max(abs(x-x_back)) < 1e-7) && (max(abs(t-t_back)) < 1e-7)
   fprintf('\tTransformed backed.\n')
else
   fprintf('\tInverse transform fail.\n')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Cosine transform
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Testing cosine transform:\n')

% Create a simple function, x = sin(pi*t) + 3*sin(4*pi*t);
% The lowest resolved frequency is *half* a cycle per unit time
x =  1 + 0*cos(pi*t/T) + 0*cos(16*(2*pi)*t/T);

[xbar, f] = CosineTransformForward(t,x);

S = T*(xbar .* conj(xbar));

df = f(2)-f(1);

% The end points only have half an increment.
x_sum = (1/T)*(sum(x(2:end-1).*x(2:end-1))*dt+x(1)*x(1)*dt/2 + x(end)*x(end)*dt/2);
S_sum = (2*S(1) + sum(S(2:end-1)) + 2*S(end))*df;

if abs(x_sum - S_sum) < 1e-7
   fprintf('\tPower matches\n')
else
   fprintf('\tPower fail\n')
end

[x_back, t_back] = CosineTransformBack(f,xbar);
if (max(abs(x-x_back)) < 1e-7) && (max(abs(t-t_back)) < 1e-7)
   fprintf('\tTransformed backed.\n')
else
   fprintf('\tInverse transform fail.\n')
end
