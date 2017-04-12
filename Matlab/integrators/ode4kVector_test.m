reps = 5000;
nDims = 1;
kappa = 1;

t = (0:100)';
deltaT = t(2)-t(1);
N = length(t);

% Notation is that f is the flux.
f = @(t,x) zeros(size(x));

x0 = zeros(reps,nDims);
x = ode4kVector(f, kappa, -Inf, Inf, t, x0);

L = 20;

figure
subplot(3,1,1)
histogram(x(:,1,1),'Normalization', 'pdf')
xlim([-L L])
ylim([0 0.3])
subplot(3,1,2)
histogram(x(:,1,25),'Normalization', 'pdf')
xlim([-L L])
ylim([0 0.3])
subplot(3,1,3)
histogram(x(:,1,50),'Normalization', 'pdf')
xlim([-L L])
ylim([0 0.3])

a = -10;
b = 10;

x = ode4kVector(f, kappa, a, Inf, t, x0);

figure
subplot(3,1,1)
histogram(x(:,1,1),'Normalization', 'pdf')
xlim([-L L])
ylim([0 0.3])
subplot(3,1,2)
histogram(x(:,1,25),'Normalization', 'pdf')
xlim([-L L])
ylim([0 0.3])
subplot(3,1,3)
histogram(x(:,1,50),'Normalization', 'pdf')
xlim([-L L])
ylim([0 0.3])

x = ode4kVector(f, kappa, -Inf, b, t, x0);

figure
subplot(3,1,1)
histogram(x(:,1,1),'Normalization', 'pdf')
xlim([-L L])
ylim([0 0.3])
subplot(3,1,2)
histogram(x(:,1,25),'Normalization', 'pdf')
xlim([-L L])
ylim([0 0.3])
subplot(3,1,3)
histogram(x(:,1,50),'Normalization', 'pdf')
xlim([-L L])
ylim([0 0.3])

x = ode4kVector(f, kappa, a, b, t, x0);

figure
subplot(3,1,1)
histogram(x(:,1,1),'Normalization', 'pdf')
xlim([-L L])
ylim([0 0.3])
subplot(3,1,2)
histogram(x(:,1,25),'Normalization', 'pdf')
xlim([-L L])
ylim([0 0.3])
subplot(3,1,3)
histogram(x(:,1,50),'Normalization', 'pdf')
xlim([-L L])
ylim([0 0.3])