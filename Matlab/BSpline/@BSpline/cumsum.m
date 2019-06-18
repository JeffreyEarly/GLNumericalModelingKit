function intspline = cumsum(spline)
%CUMSUM Indefinite integral of a BSpline


% ts = BSpline.PointsOfSupport(spline.t_knot,spline.K,0);

t_knot = cat(1,spline.t_knot(1),spline.t_knot,spline.t_knot(end));
ts = BSpline.PointsOfSupport(t_knot,spline.K+1,0);
B = BSpline.Spline(ts,t_knot,spline.K+1,1);
X = squeeze(B(:,:,1));
V = squeeze(B(:,:,2));

% This makes the first derivative match
g = spline.ValueAtPoints(ts);
V(1,:) = X(1,:);
g(1) = 0;
m = V\g;

% % Now compute the value at the LHS.
% const = X(1,:)*m;
% 
% % And subtract that value
% m1 = X\ones(size(X,1),1);
% m = m - const*m1;

intspline = BSpline(spline.K+1,t_knot,m);

