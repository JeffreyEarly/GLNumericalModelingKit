reps = 500;
nDims = 2;
kappa = 1;

t = (0:10:300)';
deltaT = t(2)-t(1);
N = length(t);

% Notation is that f is the flux.
f = @(t,x) zeros(size(x));
g = @(t,x) sqrt(2*kappa)*ones(size(x));

p0 = 0;
x0 = p0*ones(reps,nDims);

integrator = IntegratorEulerMaruyama( f, g, x0, deltaT );
pn = integrator.IntegrateAlongDimension(t);
x = squeeze(pn(:,1,:)).';
y = squeeze(pn(:,2,:)).';

D2 = x(end,:).^2 + y(end,:).^2;
kappa_out = mean(D2)/(4*t(end))

figure
plot(x,y)


% Notation is that f is the flux.
% first column, x, second column u
tau = 30;
sigma = sqrt(2*kappa)/tau;
f = @(t,x) cat(2,x(:,2),-x(:,2)/tau);
g = @(t,x) cat(2,zeros(size(x(:,1))), sigma*ones(size(x(:,2))));

x0 = zeros(reps,2);

integrator = IntegratorEulerMaruyama( f, g, x0, deltaT );
pn = integrator.IntegrateAlongDimension(t);
x = squeeze(pn(:,1,:)).';
u = squeeze(pn(:,2,:)).';

D2 = x(end,:).^2 ;
kappa_out = mean(D2)/(2*t(end))

figure
plot(t,x)