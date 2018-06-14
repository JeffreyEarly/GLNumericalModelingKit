f = @(x) sin(2*pi*x/10);

N = 64;
K = 3; % order of spline
D = K-1; % number of derivates to return
t = linspace(0,10,N)'; % observation points

for i=1:11
    N_splines = i;
    
    t_knot = BSpline.KnotPointsForSplines(t,K,N_splines);
    X = BSpline.Spline( t, t_knot, K );
    
    fprintf('(Requested, received) = (%d,%d)\n',N_splines, size(X,2));
end