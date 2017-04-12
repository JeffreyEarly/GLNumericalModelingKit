reps = 5000;
kappa = 1;

t = (0:100)';
deltaT = t(2)-t(1);
N = length(t);

% Notation is that f is the flux.
f = @(t,x) zeros(size(x));

p0 = -10;
x0 = p0*ones(reps,1);
x = ode4k(f, kappa, -Inf*ones(size(x0)), Inf*ones(size(x0)), t, x0);

L = 20;

figure
subplot(3,1,1)
histogram(x(1,:),'Normalization', 'pdf')
xlim([-L L])
ylim([0 0.3])
subplot(3,1,2)
histogram(x(25,:),'Normalization', 'pdf')
xlim([-L L])
ylim([0 0.3])
subplot(3,1,3)
histogram(x(50,:),'Normalization', 'pdf')
xlim([-L L])
ylim([0 0.3])

a = -10+p0;
b = 10+p0;

x = ode4k(f, kappa, a*ones(size(x0)), Inf*ones(size(x0)), t, x0);

figure
subplot(3,1,1)
histogram(x(1,:),'Normalization', 'pdf')
xlim([-L L])
ylim([0 0.3])
subplot(3,1,2)
histogram(x(25,:),'Normalization', 'pdf')
xlim([-L L])
ylim([0 0.3])
subplot(3,1,3)
histogram(x(50,:),'Normalization', 'pdf')
xlim([-L L])
ylim([0 0.3])

x = ode4k(f, kappa, -Inf*ones(size(x0)), b*ones(size(x0)), t, x0);

figure
subplot(3,1,1)
histogram(x(1,:),'Normalization', 'pdf')
xlim([-L L])
ylim([0 0.3])
subplot(3,1,2)
histogram(x(25,:),'Normalization', 'pdf')
xlim([-L L])
ylim([0 0.3])
subplot(3,1,3)
histogram(x(50,:),'Normalization', 'pdf')
xlim([-L L])
ylim([0 0.3])

x = ode4k(f, kappa, a*ones(size(x0)), b*ones(size(x0)), t, x0);

figure
subplot(3,1,1)
histogram(x(1,:),'Normalization', 'pdf')
xlim([-L L])
ylim([0 0.3])
subplot(3,1,2)
histogram(x(25,:),'Normalization', 'pdf')
xlim([-L L])
ylim([0 0.3])
subplot(3,1,3)
histogram(x(50,:),'Normalization', 'pdf')
xlim([-L L])
ylim([0 0.3])