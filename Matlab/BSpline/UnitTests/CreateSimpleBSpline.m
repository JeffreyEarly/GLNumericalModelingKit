K = 3;
t = (0:10)';
t_knot = BSpline.NaturalKnotsForSpline(t,K);

tq = linspace(min(t),max(t),1000)';

B = BSpline.SplineParallel( tq, t_knot, K );

figure, plot(tq,B)