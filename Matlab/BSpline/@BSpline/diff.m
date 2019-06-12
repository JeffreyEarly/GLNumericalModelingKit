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
    dspline = BSpline(spline.K,spline.t_knot,spline.m);
    dspline.K = spline.K-n;
    dspline.B = spline.B(:,:,1:(spline.K-n));
    dspline.C = spline.C(:,1:(spline.K-n));
else
    error('Can only differentiate with positive integers');
end

