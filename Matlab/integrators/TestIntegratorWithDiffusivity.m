reps = 5000;
nDims = 1;
kappa = 1;

t = (0:100)';
deltaT = t(2)-t(1);
N = length(t);

% Notation is that f is the flux.
f = @(t,x) zeros(size(x));

p0 = -10;
x0 = p0*ones(reps,nDims);
integrator = IntegratorWithDiffusivity( f, x0, deltaT, kappa, -Inf, Inf );
x = integrator.IntegrateAlongDimension(t);

Lmin = -20+p0;
Lmax = 20+p0;

figure
subplot(3,1,1)
histogram(x(:,1,1),'Normalization', 'pdf')
xlim([Lmin Lmax])
ylim([0 0.3])
subplot(3,1,2)
histogram(x(:,1,25),'Normalization', 'pdf')
xlim([Lmin Lmax])
ylim([0 0.3])
subplot(3,1,3)
histogram(x(:,1,50),'Normalization', 'pdf')
xlim([Lmin Lmax])
ylim([0 0.3])

a = -10+p0;
b = 10+p0;

integrator = IntegratorWithDiffusivity( f, x0, deltaT, kappa, a, Inf );
x = integrator.IntegrateAlongDimension(t);

figure
subplot(3,1,1)
histogram(x(:,1,1),'Normalization', 'pdf')
xlim([Lmin Lmax])
ylim([0 0.3])
subplot(3,1,2)
histogram(x(:,1,25),'Normalization', 'pdf')
xlim([Lmin Lmax])
ylim([0 0.3])
subplot(3,1,3)
histogram(x(:,1,50),'Normalization', 'pdf')
xlim([Lmin Lmax])
ylim([0 0.3])

integrator = IntegratorWithDiffusivity( f, x0, deltaT, kappa, -Inf, b );
x = integrator.IntegrateAlongDimension(t);

figure
subplot(3,1,1)
histogram(x(:,1,1),'Normalization', 'pdf')
xlim([Lmin Lmax])
ylim([0 0.3])
subplot(3,1,2)
histogram(x(:,1,25),'Normalization', 'pdf')
xlim([Lmin Lmax])
ylim([0 0.3])
subplot(3,1,3)
histogram(x(:,1,50),'Normalization', 'pdf')
xlim([Lmin Lmax])
ylim([0 0.3])

integrator = IntegratorWithDiffusivity( f, x0, deltaT, kappa, a, b );
x = integrator.IntegrateAlongDimension(t);

figure
subplot(3,1,1)
histogram(x(:,1,1),'Normalization', 'pdf')
xlim([Lmin Lmax])
ylim([0 0.3])
subplot(3,1,2)
histogram(x(:,1,25),'Normalization', 'pdf')
xlim([Lmin Lmax])
ylim([0 0.3])
subplot(3,1,3)
histogram(x(:,1,50),'Normalization', 'pdf')
xlim([Lmin Lmax])
ylim([0 0.3])