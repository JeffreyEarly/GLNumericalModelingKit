%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Fourier transform
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Testing Fourier transform:\n')
N=32; % total points
T=4.0; % total time length
t=T*(0:(N-1))'/N;

% valid frequencies (skipping zero)
omega = 2*pi*(1:floor(N/2))'/T;

% % Test all sine components (sine of nyquist is zero, so we skip it)
% for iOmega = 1:(length(omega)-1)
%     x = sin(omega(iOmega)*t);
%     Dx = DiffFourier(t,x);
%     Dx_expected = omega(iOmega)*cos(omega(iOmega)*t);
%     
%     if (max(abs(Dx-Dx_expected)) < 1e-7)
%         fprintf('\tFourier differentiation passed.\n')
%     else
%         fprintf('\tFourier differentiation failed.\n')
%     end
% end
% 
% % Test all cosine components
% for iOmega = 1:length(omega)
%     x = cos(omega(iOmega)*t);
%     Dx = DiffFourier(t,x);
%     Dx_expected = -omega(iOmega)*sin(omega(iOmega)*t);
%     
%     if (max(abs(Dx-Dx_expected)) < 1e-7)
%         fprintf('\tFourier differentiation passed.\n')
%     else
%         fprintf('\tFourier differentiation failed.\n')
%     end
% end
% 
% for iOmega = 1:(length(omega)-1)
%     x = sin(omega(iOmega)*t);
%     Dx = DiffFourier(t,x,2);
%     Dx_expected = -(omega(iOmega).^2)*sin(omega(iOmega)*t);
%     
%     if (max(abs(Dx-Dx_expected)) < 1e-7)
%         fprintf('\tFourier differentiation passed.\n')
%     else
%         fprintf('\tFourier differentiation failed.\n')
%     end
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Cosine derivates
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Testing cosine derivatives:\n')


% Test all cosine components
for iOmega = 1:length(omega)
    x = cos(omega(iOmega)*t);
    Dx = DiffCosine(t,x);
    Dx_expected = -omega(iOmega)*sin(omega(iOmega)*t);
    
    if (max(abs(Dx-Dx_expected)) < 1e-7)
        fprintf('\Cosine differentiation passed.\n')
    else
        fprintf('\Cosine differentiation failed.\n')
    end
end

return

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
   fprintf('\tPower fail\n')
end

[x_back, t_back] = CosineTransformBack(f,xbar);
if (max(abs(x-x_back)) < 1e-7) && (max(abs(t-t_back)) < 1e-7)
   fprintf('\tTransformed backed.\n')
else
   fprintf('\tInverse transform fail.\n')
end

return

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

% The sum in x is sloppy--we integrating beyond the end point, but since
% the function is zero at the end points, we end up okay.
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


