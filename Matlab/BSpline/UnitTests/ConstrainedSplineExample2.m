distribution = NormalDistribution(1);
distribution = StudentTDistribution(1,3.0);

t = linspace(0,2,11)';
f = @(x) x.^2 + 2*x + 1 + distribution.rand(size(x));
x = f(t);

K = 3;
t_knot = cat(1,min(t)*ones(K,1),max(t)*ones(K,1));

spline = ConstrainedSpline(t,x,K,t_knot,distribution,[]);
tq = linspace(min(t),max(t),10*length(t))';


figure
plot(tq,spline(tq),'LineWidth',2), hold on
scatter(t,x,6^2,'filled')