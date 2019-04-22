sigma = 1;
normal = NormalDistribution(sigma);
normal2dDistance = TwoDimDistanceDistribution(normal);

rayleigh = RayleighDistribution(sigma);

z = linspace(0,5,1000)';
figure
plot(z,normal2dDistance.pdf(z)), hold on
plot(z,rayleigh.pdf(z));

figure
subplot(2,1,1)
plot(z,1-normal2dDistance.cdf(z)./rayleigh.cdf(z)),xlog
% xlim([0 4])
subplot(2,1,2)
plot(z,1-normal2dDistance.pdf(z)./rayleigh.pdf(z)), xlog
% xlim([0 4])