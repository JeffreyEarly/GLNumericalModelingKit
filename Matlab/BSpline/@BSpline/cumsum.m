function intspline = cumsum(spline)
%CUMSUM Indefinite integral of a BSpline

m = spline.m;
K = spline.K;
t_knot = spline.t_knot;
M = length(m);

if abs(spline.x_mean) > 0
    X = spline.B(:,:,1);
    if isempty(X)
        X = BSpline.Spline( spline.t_pp, t_knot, K );
    end
    m0 = X\((spline.x_mean/spline.x_std)*ones(size(spline.t_pp)));
    m = m+m0;
end


dt = (t_knot(1+K:M+K)-t_knot(1:M))/K;
beta = zeros(length(m)+1,1);
for i=2:length(beta)
   beta(i) = beta(i-1) + m(i-1)*dt(i-1); 
end

t_knot = cat(1,spline.t_knot(1),spline.t_knot,spline.t_knot(end));
intspline = BSpline(spline.K+1,t_knot,beta);
intspline.x_std = spline.x_std;
