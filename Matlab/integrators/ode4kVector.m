function Y = ode4kVector(odefun,kappa,ymin,ymax,tspan,y0,varargin)
%ODE4  Solve differential equations with a non-adaptive method of order 4.
%   Y = ODE4(ODEFUN,TSPAN,Y0) with TSPAN = [T1, T2, T3, ... TN] integrates
%   the system of differential equations y' = f(t,y) by stepping from T0 to
%   T1 to TN. Function ODEFUN(T,Y) must return f(t,y) in a column vector.
%   The vector Y0 is the initial conditions at T0. Each row in the solution
%   array Y corresponds to a time specified in TSPAN.
%
%   Y = ODE4(ODEFUN,TSPAN,Y0,P1,P2...) passes the additional parameters
%   P1,P2... to the derivative function as ODEFUN(T,Y,P1,P2...).
%
%   This is a non-adaptive solver. The step sequence is determined by TSPAN
%   but the derivative function ODEFUN is evaluated multiple times per step.
%   The solver implements the classical Runge-Kutta method of order 4.
%
%   Example
%         tspan = 0:0.1:20;
%         y = ode4(@vdp1,tspan,[2 0]);
%         plot(tspan,y(:,1));
%     solves the system y' = vdp1(t,y) with a constant step size of 0.1,
%     and plots the first component of the solution.
%
%   This function was modified by Jeffrey J. Early to include a diffusivity
%   in any direction, as well as boundaries.

if ~isnumeric(tspan)
    error('TSPAN should be a vector of integration steps.');
end

if ~isnumeric(y0)
    error('Y0 should be a vector of initial conditions.');
end

h = diff(tspan);
if any(sign(h(1))*h <= 0)
    error('Entries of TSPAN are not in order.')
end

try
    f0 = feval(odefun,tspan(1),y0,varargin{:});
catch
    msg = ['Unable to evaluate the ODEFUN at t0,y0. ',lasterr];
    error(msg);
end

[nReps, nDims] = size(y0);

if ~isequal(size(y0),size(f0))
    error('Inconsistent sizes of Y0 and f(t0,y0).');
end

if length(kappa) == 1 && nDims > 1 % if it's a scalar, make it a vector
    kappa = kappa*ones(1,nDims);
end
if ~isequal(size(kappa,2),nDims)
    error('Inconsistent sizes of Y0 and kappa.');
end

if ~isequal([1 nDims],size(ymin)) || ~isequal([1 nDims],size(ymax))
    error('Inconsistent sizes of Y0 and ymin/ymax.');
end

f = cell(length(ymin),1);
for i = 1:length(ymin)
    if kappa(i) == 0
        % no-op
        f{i} = @(x,dt) zeros(size(x));
    elseif ymin(i) == -Inf && isfinite(ymax(i))
        % bounded above
        b = ymax(i);
        s = sqrt(2*kappa(i));
        f{i} = @(x,dt) -abs((s/sqrt(dt))*randn(size(x)) - (b-x)/dt) + (b-x)/dt;
    elseif isfinite(ymin(i)) && ymax(i) == Inf
        % bounded below
        a = ymin(i);
        s = sqrt(2*kappa(i));
        f{i} = @(x,dt) abs((s/sqrt(dt))*randn(size(x)) - (a-x)/dt) + (a-x)/dt;
    elseif isfinite(ymin(i)) && isfinite(ymax(i))
        % bounded above and below
        a = ymin(i);
        b = ymax(i);
        s = sqrt(2*kappa(i));
        f{i} = @(x,dt) abs( mod( (s/sqrt(dt))*randn(size(x))-(b-x)/dt, 2*(b-a)/dt) - (b-a)/dt ) + (a-x)/dt;
    elseif ymin(i) == -Inf && ymax(i) == Inf
        % unbounded
        f{i} = @(x,dt) sqrt(2*kappa(i)/dt)*randn(size(x));
    else
        error('Invalid diffusion boundary')
    end
end


N = length(tspan);
Y = zeros(nReps,nDims,N);
F = zeros(nReps,nDims,4);

Y(:,:,1) = y0;
for i = 2:N
    ti = tspan(i-1);
    hi = h(i-1);
    yi = Y(:,:,i-1);
    F(:,:,1) = feval(odefun,ti,yi,varargin{:});
    F(:,:,2) = feval(odefun,ti+0.5*hi,yi+0.5*hi*F(:,:,1),varargin{:});
    F(:,:,3) = feval(odefun,ti+0.5*hi,yi+0.5*hi*F(:,:,2),varargin{:});
    F(:,:,4) = feval(odefun,tspan(i),yi+hi*F(:,:,3),varargin{:});
    Y(:,:,i) = yi + (hi/6)*(F(:,:,1) + 2*F(:,:,2) + 2*F(:,:,3) + F(:,:,4));
    
    for j=1:nDims
        Y(:,j,i) = Y(:,j,i) + hi*f{j}(Y(:,j,i),hi);
    end
end
