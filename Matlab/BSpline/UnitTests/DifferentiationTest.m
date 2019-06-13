f = @(x) -x.^3 + x.^2 - 2*x + 1;
f = @(x) sin(4*pi*x);
N = 11; % number of points
K = 4; % order of spline

% first let's do a uniform grid, lower order
x = linspace(-1,1,N)';
spline = InterpolatingSpline(x,f(x),K);

% dense grid to plot the interpolant
x_dense = linspace(-1,1,10*N)';

figure
subplot(1,2,1)
plot(x_dense,f(x_dense), 'k--', 'LineWidth', 2), hold on,
plot(x_dense,spline(x_dense), 'LineWidth', 2)
scatter(x,f(x))

subplot(1,2,2)
plot(x_dense,spline(x_dense,1), 'LineWidth', 2)

% dt_knot = spline_fit.t_knot(2:end-1);
% B = BSpline.Spline(x,dt_knot,K-1);
% m = B\spline_fit(x,1);
% dspline = BSpline(K-1,dt_knot,m);

dspline = diff(spline,1);

hold on
plot(x_dense,dspline(x_dense), 'LineWidth', 2)


figure
plot(x_dense,spline(x_dense,3), 'LineWidth', 2), hold on
plot(x_dense,dspline(x_dense,2)), hold on