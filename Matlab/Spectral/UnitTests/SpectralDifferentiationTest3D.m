
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Fourier derivatives
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Testing Fourier transform:\n')

Nx = 4;
Ny = 8;
Nz = 16;

Lx = 1.0;
Ly = 10.0;
Lz = 4.0;

x = (Lx/Nx)*(0:Nx-1)';
y = (Ly/Ny)*(0:Ny-1)';
z = (Lz/Nz)*(0:Nz-1)';

omega_x = 2*pi/Lx;

% valid frequencies (skipping zero)
omega_z = 2*pi*(1:floor(Nz/2))'/Lz;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1 derivative, cosine, 1st-dimension
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Z,X,Y] = ndgrid(z,x,y);

f = @(omega) cos(omega_x*X).*cos(omega*Z);
Df_analytical = @(omega) -omega*cos(omega_x*X).*sin(omega*Z);
Df_numerical = @(x) DiffFourier(z,x);
testname = 'Fourier differentiation of cosine (numDerivs=1), 1st dimension ';

ReportErrors(f,Df_analytical,Df_numerical,testname,omega_z);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1 derivative, cosine, 2nd-dimension
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[X,Z,Y] = ndgrid(x,z,y);

f = @(omega) cos(omega_x*X).*cos(omega*Z);
Df_analytical = @(omega) -omega*cos(omega_x*X).*sin(omega*Z);
Df_numerical = @(x) DiffFourier(z,x);
testname = 'Fourier differentiation of cosine (numDerivs=1), 2nd dimension ';

ReportErrors(f,Df_analytical,Df_numerical,testname,omega_z);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1 derivative, cosine, 3rd-dimension
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[X,Y,Z] = ndgrid(x,y,z);

f = @(omega) cos(omega_x*X).*cos(omega*Z);
Df_analytical = @(omega) -omega*cos(omega_x*X).*sin(omega*Z);
Df_numerical = @(x) DiffFourier(z,x);
testname = 'Fourier differentiation of cosine (numDerivs=1), 3rd dimension ';

ReportErrors(f,Df_analytical,Df_numerical,testname,omega_z);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Cosine derivatives
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Testing cosine derivatives:\n')

Nx = 4;
Ny = 8;
Nz = 16;

Lx = 1.0;
Ly = 10.0;
Lz = 4.0;

x = (Lx/Nx)*(0:Nx-1)';
y = (Ly/Ny)*(0:Ny-1)';
z = Lz*(0:(Nz-1))'/(Nz-1);

omega_x = 2*pi/Lx;

% valid frequencies (skipping zero)
df = 1/((Nz-1)*(z(2)-z(1)));
omega_z = pi*df*(0:(Nz-1))';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1 derivative, cosine, 1st-dimension
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Z,X,Y] = ndgrid(z,x,y);

f = @(omega) cos(omega_x*X).*cos(omega*Z);
Df_analytical = @(omega) -omega*cos(omega_x*X).*sin(omega*Z);
Df_numerical = @(x) DiffCosine(z,x);
testname = 'Cosine differentiation (numDerivs=1), 1st dimension ';

ReportErrors(f,Df_analytical,Df_numerical,testname,omega_z);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1 derivative, cosine, 2nd-dimension
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[X,Z,Y] = ndgrid(x,z,y);

f = @(omega) cos(omega_x*X).*cos(omega*Z);
Df_analytical = @(omega) -omega*cos(omega_x*X).*sin(omega*Z);
Df_numerical = @(x) DiffCosine(z,x);
testname = 'Cosine differentiation (numDerivs=1), 2nd dimension ';

ReportErrors(f,Df_analytical,Df_numerical,testname,omega_z);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1 derivative, cosine, 3rd-dimension
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[X,Y,Z] = ndgrid(x,y,z);

f = @(omega) cos(omega_x*X).*cos(omega*Z);
Df_analytical = @(omega) -omega*cos(omega_x*X).*sin(omega*Z);
Df_numerical = @(x) DiffCosine(z,x);
testname = 'Cosine differentiation (numDerivs=1), 3rd dimension ';

ReportErrors(f,Df_analytical,Df_numerical,testname,omega_z);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Sine derivatives
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Testing sine derivatives:\n')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1 derivative, sine, 1st-dimension
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Z,X,Y] = ndgrid(z,x,y);

f = @(omega) cos(omega_x*X).*sin(omega*Z);
Df_analytical = @(omega) omega*cos(omega_x*X).*cos(omega*Z);
Df_numerical = @(x) DiffSine(z,x);
testname = 'Sine differentiation (numDerivs=1), 1st dimension ';

ReportErrors(f,Df_analytical,Df_numerical,testname,omega_z(1:end-1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1 derivative, sine, 2nd-dimension
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[X,Z,Y] = ndgrid(x,z,y);

f = @(omega) cos(omega_x*X).*sin(omega*Z);
Df_analytical = @(omega) omega*cos(omega_x*X).*cos(omega*Z);
Df_numerical = @(x) DiffSine(z,x);
testname = 'Sine differentiation (numDerivs=1), 2nd dimension ';

ReportErrors(f,Df_analytical,Df_numerical,testname,omega_z(1:end-1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1 derivative, sine, 3rd-dimension
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[X,Y,Z] = ndgrid(x,y,z);

f = @(omega) cos(omega_x*X).*sin(omega*Z);
Df_analytical = @(omega) omega*cos(omega_x*X).*cos(omega*Z);
Df_numerical = @(x) DiffSine(z,x);
testname = 'Sine differentiation (numDerivs=1), 3rd dimension ';

ReportErrors(f,Df_analytical,Df_numerical,testname,omega_z(1:end-1));

return

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

