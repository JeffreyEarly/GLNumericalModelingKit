function dspline = diff(spline,n)
%DIFF Differentiation of a BSpline

if nargin < 2
    n = 1;
end

if n == 0
    dspline = spline;
elseif n >= spline.K
    dspline = BSpline(spline.K,spline.t_knot,zeros(size(spline.m)));
    n = spline.K-1;
    dspline.K = spline.K-n;
    dspline.B = spline.B(:,:,1:(spline.K-n));
    dspline.C = zeros(size(spline.C(:,1:(spline.K-n))));
elseif ( ( n > 0 ) && ( round(n) == n ) )   % Positive integer
    %     dspline = BSpline(spline.K-1,spline.t_knot(2:end-1),spline.m(1:end-1));
    %     dspline.K = spline.K-n;
    %     dspline.B = spline.B(:,:,(n+1):(spline.K-n));
    %     dspline.C = spline.C(:,1:(spline.K-n));
    
    ts = BSpline.PointsOfSupport(spline.t_knot,spline.K,0);
    t_knot = spline.t_knot(2:end-1);
    X = BSpline.Spline(ts,t_knot,spline.K-1);
    
    g = spline.ValueAtPoints(ts,1);
    m = X\g;
    dspline = BSpline(spline.K-1,t_knot,m);
else
    error('Can only differentiate with positive integers');
end

