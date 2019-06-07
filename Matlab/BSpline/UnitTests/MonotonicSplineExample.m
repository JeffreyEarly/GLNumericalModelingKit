rng(1)
distribution = NormalDistribution(0.2);
% distribution = StudentTDistribution(1,3.0);

t = linspace(0,2,11)';
f = @(x) tanh(x-1) + distribution.rand(size(x));
x = f(t);

K = 3;
t_knot = cat(1,min(t)*ones(K,1),[0.5; 1; 1.5],max(t)*ones(K,1));

spline0 = ConstrainedSpline(t,x,K,t_knot,distribution,[]);
spline = ConstrainedMonotonicSpline(t,x,K,t_knot,distribution,[]);
tq = linspace(min(t),max(t),10*length(t))';


figure
plot(tq,spline(tq),'LineWidth',2), hold on
plot(tq,spline0(tq),'LineWidth',2)
scatter(t,x,6^2,'filled')
legend('monotonic','unconstrained')