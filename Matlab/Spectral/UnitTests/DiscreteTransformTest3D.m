

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Testing cosine derivatives:\n')

Nx = 8;
Ny = 32;
Nz = 16+1;

Lx = 1.0;
Ly = 10.0;
Lz = 4.0;

x = (Lx/Nx)*(0:Nx-1)';
y = (Ly/Ny)*(0:Ny-1)';
z = Lz*(0:(Nz-1))'/(Nz-1);

[X,Y,Z] = ndgrid(x,y,z);

omega = 2*pi*(0:floor(Nx/2))'/Lx;

f = @(omega) cos(omega*X).*cos(pi*Z);

u = f(omega(4));

[ubar,k] = FourierTransformForward(x,u);
[ubar,m] = CosineTransformForward(x,ubar);
u_back = FourierTransformBack(

return

Df_analytical = @(omega) -omega*sin(omega*X).*cos(pi*Z);
Df_numerical = @(val) DiffFourier(x,val);
testname = 'Fourier differentiation of cosine (numDerivs=1), 1st dimension';

ReportErrors(f,Df_analytical,Df_numerical,testname,omega);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Cosine derivatives
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Testing cosine derivatives:\n')

Nx = 4;
Ny = 8;
Nz = 16+1;

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
testname = 'Cosine differentiation (numDerivs=1), 1st dimension';

ReportErrors(f,Df_analytical,Df_numerical,testname,omega_z);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1 derivative, cosine, 2nd-dimension
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[X,Z,Y] = ndgrid(x,z,y);

f = @(omega) cos(omega_x*X).*cos(omega*Z);
Df_analytical = @(omega) -omega*cos(omega_x*X).*sin(omega*Z);
Df_numerical = @(x) DiffCosine(z,x);
testname = 'Cosine differentiation (numDerivs=1), 2nd dimension';

ReportErrors(f,Df_analytical,Df_numerical,testname,omega_z);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1 derivative, cosine, 3rd-dimension
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[X,Y,Z] = ndgrid(x,y,z);

f = @(omega) cos(omega_x*X).*cos(omega*Z);
Df_analytical = @(omega) -omega*cos(omega_x*X).*sin(omega*Z);
Df_numerical = @(x) DiffCosine(z,x);
testname = 'Cosine differentiation (numDerivs=1), 3rd dimension';

ReportErrors(f,Df_analytical,Df_numerical,testname,omega_z);

