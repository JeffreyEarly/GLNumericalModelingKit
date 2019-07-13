function intspline = cumsum(spline)
%CUMSUM Indefinite integral of a BSpline

m = spline.m;
K = spline.K;
t_knot = spline.t_knot;
M = length(m);

if abs(spline.x_mean) > 0 || abs(spline.x_std - 1) > 0
%     X = spline.B(:,:,1);
%     if isempty(X)
%         X = BSpline.Spline( spline.t_pp, t_knot, K );
%     end
    t = BSpline.PointsOfSupport(spline.t_knot,spline.K);
    X = BSpline.Spline(t,spline.t_knot,spline.K);
    m = spline.x_std*spline.m + X\(spline.x_mean*ones(length(t),1));
end


dt = (t_knot(1+K:M+K)-t_knot(1:M))/K;
beta = zeros(length(m)+1,1);
for i=2:length(beta)
   beta(i) = beta(i-1) + m(i-1)*dt(i-1); 
end

t_knot = cat(1,spline.t_knot(1),spline.t_knot,spline.t_knot(end));
intspline = BSpline(spline.K+1,t_knot,beta);
% intspline.x_std = spline.x_std;
