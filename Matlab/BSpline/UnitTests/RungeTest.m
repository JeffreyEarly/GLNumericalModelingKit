f = @(x) 1./(25*x.*x+1);
N = 11; % number of points
K = 4; % order of spline

x = linspace(-1,1,N)';

spline_fit = BSpline(x,f(x),K);

x_dense = linspace(-1,1,10*N)';

figure
plot(x_dense,f(x_dense), 'k--', 'LineWidth', 2), hold on,
plot(x_dense,spline_fit(x_dense))


% x = [-1 -0.33 0.33 1]';
% y = [0 1 2 0]';
% 
% spline_fit = BSpline(x,y,4);
% interpolant = griddedInterpolant(x,y,'spline');
% figure
% scatter(x,y), hold on,
% plot(x_dense,spline_fit(x_dense))
% plot(x_dense,interpolant(x_dense))
