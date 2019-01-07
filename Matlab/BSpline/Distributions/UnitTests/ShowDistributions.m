pct = 0.05;
distribution = AddedDistribution(pct,NormalDistribution(800),StudentTDistribution(8.5,4.5));

z = linspace(-500,500,1000);
% figure, plot(z,distribution.pdf(z))

distDistribution = TwoDimDistanceDistribution(distribution);
figure, plot(z,distDistribution.cdf(z))