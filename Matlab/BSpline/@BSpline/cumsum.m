function intspline = cumsum(spline)
%CUMSUM Indefinite integral of a BSpline

ts = BSpline.PointsOfSupport(spline.t_knot,spline.K,0);
t_knot = cat(1,spline.t_knot(1),spline.t_knot,spline.t_knot(end));
B = BSpline.Spline(ts,t_knot,spline.K+1,1);
V = squeeze(B(:,:,2));

g = spline.ValueAtPoints(ts);
m = V\g;
intspline = BSpline(spline.K+1,t_knot,m);

