sigma = 1;
normal = NormalDistribution(sigma);
normal2dDistance = TwoDimDistanceDistribution(normal);

rayleigh = RayleighDistribution(sigma);

z = linspace(0,10,1000)';
figure
plot(z,normal2dDistance.pdf(z)), hold on
plot(z,rayleigh.pdf(z));

figure
plot(z,normal2dDistance.cdf(z)-rayleigh.cdf(z))
xlim([0 4])