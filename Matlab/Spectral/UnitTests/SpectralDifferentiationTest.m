%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Fourier derivatives
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Testing Fourier transform:\n')
N=32; % total points
T=4.0; % total time length
t=T*(0:(N-1))'/N;

% valid frequencies (skipping zero)
omega = 2*pi*(1:floor(N/2))'/T;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1 derivative, sine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f = @(omega) sin(omega*t);
Df_analytical = @(omega) omega*cos(omega*t);
Df_numerical = @(x) DiffFourier(t,x);
testname = 'Fourier differentiation of sine (numDerivs=1) ';

ReportErrors(f,Df_analytical,Df_numerical,testname,omega(1:end-1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1 derivative, cosine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f = @(omega) cos(omega*t);
Df_analytical = @(omega) -omega*sin(omega*t);
Df_numerical = @(x) DiffFourier(t,x);
testname = 'Fourier differentiation of cosine (numDerivs=1) ';

ReportErrors(f,Df_analytical,Df_numerical,testname,omega);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2 derivatives, sine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f = @(omega) sin(omega*t);
Df_analytical = @(omega) -(omega^2)*sin(omega*t);
Df_numerical = @(x) DiffFourier(t,x,2);
testname = 'Fourier differentiation of sine (numDerivs=2) ';

ReportErrors(f,Df_analytical,Df_numerical,testname,omega(1:end-1));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Cosine derivatives
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Testing cosine derivatives:\n')

N=33; % total points
T=4.0; % total time length
t=T*(0:(N-1))'/(N-1);

% valid frequencies (skipping zero)
df = 1/((N-1)*(t(2)-t(1)));
omega = pi*df*(0:(N-1))';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1 derivative
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f = @(omega) cos(omega*t);
Df_analytical = @(omega) -omega*sin(omega*t);
Df_numerical = @(x) DiffCosine(t,x);
testname = 'Cosine differentiation (numDerivs=1) ';

ReportErrors(f,Df_analytical,Df_numerical,testname,omega);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2 derivatives
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f = @(omega) cos(omega*t);
Df_analytical = @(omega) -(omega^2)*cos(omega*t);
Df_numerical = @(x) DiffCosine(t,x,2);
testname = 'Cosine differentiation (numDerivs=2) ';

ReportErrors(f,Df_analytical,Df_numerical,testname,omega);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3 derivatives
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f = @(omega) cos(omega*t);
Df_analytical = @(omega) (omega^3)*sin(omega*t);
Df_numerical = @(x) DiffCosine(t,x,3);
testname = 'Cosine differentiation (numDerivs=3) ';

ReportErrors(f,Df_analytical,Df_numerical,testname,omega);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4 derivatives
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f = @(omega) cos(omega*t);
Df_analytical = @(omega) (omega^4)*cos(omega*t);
Df_numerical = @(x) DiffCosine(t,x,4);
testname = 'Cosine differentiation (numDerivs=4) ';

ReportErrors(f,Df_analytical,Df_numerical,testname,omega);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Sine derivatives
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Testing sine derivatives:\n')

N=33; % total points
T=4.0; % total time length
t=T*(0:(N-1))'/(N-1);

% valid frequencies (skipping zero)
df = 1/((N-1)*(t(2)-t(1)));
omega = pi*df*(0:(N-1))';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1 derivative
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f = @(omega) sin(omega*t);
Df_analytical = @(omega) omega*cos(omega*t);
Df_numerical = @(x) DiffSine(t,x);
testname = 'Sine differentiation (numDerivs=1) ';

ReportErrors(f,Df_analytical,Df_numerical,testname,omega(1:end-1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2 derivatives
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f = @(omega) sin(omega*t);
Df_analytical = @(omega) -(omega^2)*sin(omega*t);
Df_numerical = @(x) DiffSine(t,x,2);
testname = 'Sine differentiation (numDerivs=2) ';

ReportErrors(f,Df_analytical,Df_numerical,testname,omega);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3 derivatives
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f = @(omega) sin(omega*t);
Df_analytical = @(omega) -(omega^3)*cos(omega*t);
Df_numerical = @(x) DiffSine(t,x,3);
testname = 'Sine differentiation (numDerivs=3) ';

ReportErrors(f,Df_analytical,Df_numerical,testname,omega(1:end-1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4 derivatives
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f = @(omega) sin(omega*t);
Df_analytical = @(omega) (omega^4)*sin(omega*t);
Df_numerical = @(x) DiffSine(t,x,4);
testname = 'Sine differentiation (numDerivs=4) ';

ReportErrors(f,Df_analytical,Df_numerical,testname,omega);

