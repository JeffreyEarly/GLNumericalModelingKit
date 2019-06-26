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
    D = n;
    m = spline.m;
    K = spline.K;
    t_knot = spline.t_knot;
    M = length(m);
    
    alpha = zeros(length(m),D+1);
    alpha(:,1) = m; % first column is the existing coefficients
    
    for d=1:D
        dm = diff(alpha(:,d));
        dt = (t_knot(1+K-d:M+K-d)-t_knot(1:M))/(K-d);
        alpha(1:end-d,d+1) = dm./dt(d+1:end);
    end
    
    dspline = BSpline(K-D,t_knot((1+D):(end-D)),alpha(1:end-D,D+1));
    dspline.x_std = spline.x_std;   
else
    error('Can only differentiate with positive integers');
end

