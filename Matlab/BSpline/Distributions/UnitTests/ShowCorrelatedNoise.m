sigma = 2;
tau = 10;
distribution = NormalDistribution(sigma);
distribution.rho = @(z) exp(-(z/tau).^2);

t = (0:500).';
z = distribution.noise(t);

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