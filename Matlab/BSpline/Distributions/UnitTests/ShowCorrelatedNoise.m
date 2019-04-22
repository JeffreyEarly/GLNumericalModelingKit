%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Create a figure with uncorrelated noise
%
sigma = 1;
distribution = NormalDistribution(sigma);

n = 1000;
z = distribution.rand([n 1]);

figure
histogram(z,50,'Normalization','pdf')
zdist = linspace(min(z),max(z),100)';
hold on
plot(zdist,distribution.pdf(zdist),'LineWidth', 2)

print('-depsc2', '../figures/normaldistribution.eps')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Create a figure with *correlated* noise
%

tau = 10;
distribution.rho = @(z) exp(-(z/tau).^2);

t = (0:500).';
z = distribution.noise(t);

figure
plot(t,z)

print('-depsc2', '../figures/correlatednoise.eps')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Create a figure with a bunch of different metrics
%

% compute the autocorrelation
tau = t(1:find(t<3*tau,1,'last'));
[ac,dof] = Autocorrelation(z,length(tau)-1);

% estimate the standard error
SE_indep = tau(2:end);
SE =  sqrt((1 + 2*cumsum(ac.^2))./dof);
SE(1) = sqrt(1/dof(1)); % first point is lag 1
SE(end) = []; % there is no end point

figure

subplot(2,2,1)
plot(t,z)
title(sprintf('mean=%.2f, std=%.2f',mean(z),std(z)));
xlabel('time (s)')

subplot(2,2,3)
plot(tau,ac), hold on
plot(tau,distribution.rho(tau))
plot(SE_indep, [3*SE,-3*SE], 'LineWidth', 1.5, 'Color',0.4*[1.0 1.0 1.0] )
legend('actual','theoretical','3\sigma standard error')
xlabel('lag (s)')
ylabel('autocorrelation')
ylim([0 1])

subplot(2,2,2)
histogram(z,100,'Normalization', 'pdf'), hold on
zdist = linspace(min(z),max(z),100)';
plot(zdist,distribution.pdf(zdist))

subplot(2,2,4)
histogram(z,100,'Normalization', 'cdf'), hold on
plot(zdist,distribution.cdf(zdist))