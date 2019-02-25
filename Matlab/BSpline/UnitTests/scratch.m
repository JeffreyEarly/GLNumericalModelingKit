L = 10;
f = @(x) sin(2*pi*x/L);

N = 10;
x = linspace(0,2*L,N)';

sigma = 0.1;
y = f(x);

x_dense = linspace(0,2*L,10*N)';

figure
plot(x_dense,f(x_dense)), hold on
scatter(x,y,'k')
legend('true function', 'observations')

spline = TensionSpline(x,y,NormalDistribution(sigma),'S',3,'T',3);

figure
plot(x_dense,f(x_dense)), hold on
scatter(x,y,'k')
plot(x_dense,spline(x_dense),'LineWidth',2)
legend('true function', 'noisy data', 'tension spline fit')