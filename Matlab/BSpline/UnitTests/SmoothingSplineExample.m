sigma = 0.25;
distribution = NormalDistribution(sigma);
% distribution = StudentTDistribution(1,3.0);


z = linspace(-4000,0,31)';
N0 = 5.2e-3; % reference buoyancy frequency, radians/seconds
g = 9.81;
rho0 = 1025;
L_gm = 1300;
rhoFunc = @(z) rho0*(1 + L_gm*N0*N0/(2*g)*(1 - exp(2*z/L_gm))) + distribution.rand(size(z));

rho = rhoFunc(z);

spline0 = SmoothingSpline(z,rho,distribution);
spline = MonotonicSmoothingSpline(z,rho,distribution);


% K = 3;
% t_knot = cat(1,min(z)*ones(K,1),max(z)*ones(K,1));
% spline = ConstrainedMonotonicSpline(z,rho,K,t_knot,distribution,[]);

zq = linspace(min(z),max(z),10*length(z))';

figure
subplot(1,2,1)
plot(spline0(zq),zq,'LineWidth',2), hold on
plot(spline(zq),zq,'LineWidth',2)
scatter(rho,z,6^2,'filled')

subplot(1,2,2)
plot(-(g/rho0)*spline0(zq,1),zq,'LineWidth',2), hold on
plot(-(g/rho0)*spline(zq,1),zq,'LineWidth',2)