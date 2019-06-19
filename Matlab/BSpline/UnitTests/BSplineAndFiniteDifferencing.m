% I find that a K=3 spline looks like 5th order accuracy finite
% differencing in the center of the domain, and 2.5th order near the
% boundary.
K = 3;
t = linspace(0,20,21)';
t_knot = InterpolatingSpline.KnotPointsForPoints( t, K );

% t_knot(end-2*K-1:end-K-1) = [];
% t_knot(K+1:2*K+1) = [];

B = BSpline.Spline( t, t_knot, K, 2 );
X = B(:,:,1);
V = B(:,:,2);
A = B(:,:,3);

iX = inv(X);
D = A*iX;

f = @(x) x.^2;
m = X\f(t);
df = V*m;

