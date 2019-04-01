% create some 'true' function
L = 10;
f = @(x) sin(2*pi*x/L);

% define a grid of observation points
N = 30;
x = linspace(0,L,N)';

% now sample the function at the observation points
sigma = 0.1;
y = f(x) + sigma*randn(N,1);

% create a dense grid of points to visualize the function 
x_dense = linspace(0,L,10*N)';

figure
plot(x_dense,f(x_dense)), hold on
scatter(x,y,'k')
legend('true function', 'noisy data')
% print('-depsc2', '../figures/noisydata.eps')

% now create a tension spline fit to the data
spline = TensionSpline(x,y,NormalDistribution(sigma));

figure
plot(x_dense,f(x_dense)), hold on
scatter(x,y,'k')
plot(x_dense,spline(x_dense),'LineWidth',2)
legend('true function', 'noisy data', 'tension spline fit')
% print('-depsc2', '../figures/noisydatawithtensionspline.eps')