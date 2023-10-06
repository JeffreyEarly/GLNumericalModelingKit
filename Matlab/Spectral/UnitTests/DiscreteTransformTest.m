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

[xbar, f] = FourierTransformForward(t,x);

S = T*(xbar .* conj(xbar));

dt = t(2)-t(1);
df = f(2)-f(1);

x_sum = (1/T)*sum(x.*x)*dt;
S_sum = sum(S)*df;

if abs(x_sum - S_sum) < 1e-7
   fprintf('\tPower matches.\n')
else
   fprintf('\tPower FAIL\n')
end

[x_back, t_back] = FourierTransformBack(f,xbar,1);
if (max(abs(x-x_back)) < 1e-7) && (max(abs(t-t_back)) < 1e-7)
   fprintf('\tTransformed backed.\n')
else
   fprintf('\tInverse transform FAIL.\n')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Sine transform
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Testing real sine transform:\n')

N=33; % total points
T=4.0; % total time length
t=T*(0:(N-1))'/(N-1);

% Create a simple function, x = sin(pi*t) + 3*sin(4*pi*t);
% The lowest resolved frequency is *half* a cycle per unit time
x =  sin(pi*t/T) + 2*sin(8*(2*pi)*t/T) + 3*sin(15.5*(2*pi)*t/T);

[xbar,f] = SineTransformForward(t,x);

S = T*(xbar .* conj(xbar));

df = f(2)-f(1);

% The sum in x is sloppy--we integrating beyond the end point, but since
% the function is zero at the end points, we end up okay.
x_sum = (1/T)*sum(x.*x)*dt;
S_sum = sum(S)*df;

if abs(x_sum - S_sum) < 1e-7
   fprintf('\tPower matches\n')
else
   fprintf('\tPower FAIL\n')
end

[x_back, t_back] = SineTransformBack(f,xbar);
if (max(abs(x-x_back)) < 1e-7) && (max(abs(t-t_back)) < 1e-7)
   fprintf('\tTransformed backed.\n')
else
   fprintf('\tInverse transform FAIL.\n')
end

if (~isreal(xbar))
   fprintf('\tTransform forward is not purely real\n'); 
end

if (~isreal(x_back))
   fprintf('\tTransform back is not purely real\n'); 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Sine transform (Imaginary)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Testing imaginary sine transform:\n')

N=33; % total points
T=4.0; % total time length
t=T*(0:(N-1))'/(N-1);

% Create a simple function, x = sin(pi*t) + 3*sin(4*pi*t);
% The lowest resolved frequency is *half* a cycle per unit time
x =  sqrt(-1)*(sin(pi*t/T) + 2*sin(8*(2*pi)*t/T) + 3*sin(15.5*(2*pi)*t/T));

[xbar,f] = SineTransformForward(t,x);

S = T*(xbar .* conj(xbar));

df = f(2)-f(1);

% The sum in x is sloppy--we integrating beyond the end point, but since
% the function is zero at the end points, we end up okay.
x_sum = (1/T)*sum(x.*conj(x))*dt;
S_sum = sum(S)*df;

if abs(x_sum - S_sum) < 1e-7
   fprintf('\tPower matches\n')
else
   fprintf('\tPower FAIL\n')
end

[x_back, t_back] = SineTransformBack(f,xbar);
if (max(abs(x-x_back)) < 1e-7) && (max(abs(t-t_back)) < 1e-7)
   fprintf('\tTransformed backed.\n')
else
   fprintf('\tInverse transform FAIL.\n')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Sine transform (Complex)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Testing complex sine transform:\n')

N=33; % total points
T=4.0; % total time length
t=T*(0:(N-1))'/(N-1);

% Create a simple function, x = sin(pi*t) + 3*sin(4*pi*t);
% The lowest resolved frequency is *half* a cycle per unit time
x =  (1+sqrt(-1))*sin(pi*t/T) + (2-3*sqrt(-1))*sin(8*(2*pi)*t/T) + (3+sqrt(-1))*sin(15.5*(2*pi)*t/T);

[xbar,f] = SineTransformForward(t,x);

S = T*(xbar .* conj(xbar));

df = f(2)-f(1);

% The sum in x is sloppy--we integrating beyond the end point, but since
% the function is zero at the end points, we end up okay.
x_sum = (1/T)*sum(x.*conj(x))*dt;
S_sum = sum(S)*df;

if abs(x_sum - S_sum) < 1e-7
   fprintf('\tPower matches\n')
else
   fprintf('\tPower FAIL\n')
end

[x_back, t_back] = SineTransformBack(f,xbar);
if (max(abs(x-x_back)) < 1e-7) && (max(abs(t-t_back)) < 1e-7)
   fprintf('\tTransformed backed.\n')
else
   fprintf('\tInverse transform FAIL.\n')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Cosine transform (Real)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Testing real cosine transform:\n')

% Create a simple function, x = sin(pi*t) + 3*sin(4*pi*t);
% The lowest resolved frequency is *half* a cycle per unit time
x =  1 + 2*cos(pi*t/T) + 3*cos(16*(2*pi)*t/T);

[xbar, f] = CosineTransformForward(t,x);

S = T*(xbar .* conj(xbar));

df = f(2)-f(1);

% The end points only have half an increment.
x_sum = (1/T)*(sum(x(2:end-1).*x(2:end-1))*dt+x(1)*x(1)*dt/2 + x(end)*x(end)*dt/2);
S_sum = ( S(1)/2 + sum(S(2:end-1)) + 2*S(end))*df;

if abs(x_sum - S_sum) < 1e-7
   fprintf('\tPower matches\n')
else
   fprintf('\tPower FAIL\n')
end

[x_back, t_back] = CosineTransformBack(f,xbar);
if (max(abs(x-x_back)) < 1e-7) && (max(abs(t-t_back)) < 1e-7)
   fprintf('\tTransformed backed.\n')
else
   fprintf('\tInverse transform FAIL.\n')
end

if (~isreal(xbar))
   fprintf('\tTransform forward is not purely real\n'); 
end

if (~isreal(x_back))
   fprintf('\tTransform back is not purely real\n'); 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Cosine transform (Imaginary)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Testing imaginary cosine transform:\n')

% Create a simple function, x = sin(pi*t) + 3*sin(4*pi*t);
% The lowest resolved frequency is *half* a cycle per unit time
x =  sqrt(-1)*(1 + 2*cos(pi*t/T) + 3*cos(16*(2*pi)*t/T));

[xbar, f] = CosineTransformForward(t,x);

S = T*(xbar .* conj(xbar));

df = f(2)-f(1);

% The end points only have half an increment.
x_sum = (1/T)*(sum(x(2:end-1).*conj(x(2:end-1)))*dt+x(1)*conj(x(1))*dt/2 + x(end)*conj(x(end))*dt/2);
S_sum = ( S(1)/2 + sum(S(2:end-1)) + 2*S(end))*df;

if abs(x_sum - S_sum) < 1e-7
   fprintf('\tPower matches\n')
else
   fprintf('\tPower FAIL\n')
end

[x_back, t_back] = CosineTransformBack(f,xbar);
if (max(abs(x-x_back)) < 1e-7) && (max(abs(t-t_back)) < 1e-7)
   fprintf('\tTransformed backed.\n')
else
   fprintf('\tInverse transform FAIL.\n')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Cosine transform (Complex)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Testing Complex cosine transform:\n')

% Create a simple function, x = sin(pi*t) + 3*sin(4*pi*t);
% The lowest resolved frequency is *half* a cycle per unit time
x =  (1+2*sqrt(-1) + (2+3*sqrt(-1))*cos(pi*t/T) + (3-1*sqrt(-1))*cos(16*(2*pi)*t/T));

[xbar, f] = CosineTransformForward(t,x);

S = T*(xbar .* conj(xbar));

df = f(2)-f(1);

% The end points only have half an increment.
x_sum = (1/T)*(sum(x(2:end-1).*conj(x(2:end-1)))*dt+x(1)*conj(x(1))*dt/2 + x(end)*conj(x(end))*dt/2);
S_sum = ( S(1)/2 + sum(S(2:end-1)) + 2*S(end))*df;

if abs(x_sum - S_sum) < 1e-7
   fprintf('\tPower matches\n')
else
   fprintf('\tPower FAIL\n')
end

[x_back, t_back] = CosineTransformBack(f,xbar);
if (max(abs(x-x_back)) < 1e-7) && (max(abs(t-t_back)) < 1e-7)
   fprintf('\tTransformed backed.\n')
else
   fprintf('\tInverse transform FAIL.\n')
end

