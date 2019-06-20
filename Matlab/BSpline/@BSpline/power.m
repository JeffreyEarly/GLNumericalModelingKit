function splineb = power(spline,b)
%POWER Power of a BSpline

if b == 1
    splineb = spline;
else
    ts = BSpline.PointsOfSupport(spline.t_knot,spline.K,0);
    g = spline.ValueAtPoints(ts);
    g(abs(2*eps)>g) = 0;
%     splineb = InterpolatingSpline(ts,(g.^b),ceil(b*spline.K));
    K = ceil(b*spline.K);
    t_knot = InterpolatingSpline.KnotPointsForPoints(ts,K);
    splineb = ConstrainedSpline(ts,(g.^b),K,t_knot,[],struct('global',ShapeConstraint.positive));
%     X = BSpline.Spline(ts,spline.t_knot,spline.K);
%     g = spline.ValueAtPoints(ts);
%     m = X\(g.^b);
%     splineb = BSpline(spline.K,spline.t_knot,m);
end

