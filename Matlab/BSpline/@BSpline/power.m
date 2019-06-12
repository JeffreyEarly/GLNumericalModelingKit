function splineb = power(spline,b)
%POWER Power of a BSpline

if b == 1
    splineb = spline;
else
    ts = BSpline.PointsOfSupport(spline.t_knot,spline.K,0);
    X = BSpline.Spline(ts,spline.t_knot,spline.K);
    g = spline.ValueAtPoints(ts);
    m = X\(g.^b);
    splineb = BSpline(spline.K,spline.t_knot,m);
end

