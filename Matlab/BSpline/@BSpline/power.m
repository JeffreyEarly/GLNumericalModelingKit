function splineb = power(spline,b,constraints)
%POWER Power of a BSpline

if b == 1
    splineb = spline;
else
    ts = BSpline.PointsOfSupport(spline.t_knot,spline.K,0);
    g = spline.ValueAtPoints(ts);
    g(abs(2*eps)>g) = 0;
%         splineb = InterpolatingSpline(ts,(g.^b),ceil(b*spline.K));
        K = ceil(b*spline.K);
%     K = spline.K;
    t_knot = InterpolatingSpline.KnotPointsForPoints(ts,K);
    if exist('constraints','var') && ~isempty(constraints)
        splineb = ConstrainedSpline(ts,(g.^b),K,t_knot,[],constraints);
    else
        X = BSpline.Spline(ts,t_knot,K);
        m = X\(g.^b);
        splineb = BSpline(K,t_knot,m);
    end
    
end

