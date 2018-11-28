f = @(x) 1./(25*x.*x+1);
N = 11; % number of points
K = 4; % order of spline

% first let's do a uniform grid, lower order
x = linspace(-1,1,N)';
spline_fit = InterpolatingSpline(x,f(x),K);

% dense grid to plot the interpolant
x_dense = linspace(-1,1,10*N)';

figure
subplot(1,2,1)
plot(x_dense,f(x_dense), 'k--', 'LineWidth', 2), hold on,
plot(x_dense,spline_fit(x_dense), 'LineWidth', 2)
scatter(x,f(x))
maxError = max( abs( f(x_dense)-spline_fit(x_dense) ) );
title(sprintf('uniform grid, N=%d, K=%d, max error = %.2g',N,K,maxError))

% first let's do a chebyshev grid, high order
K = 4;
x = cos((0:N-1)'*pi/(N-1));
spline_fit = InterpolatingSpline(x,f(x),K);

subplot(1,2,2)
plot(x_dense,f(x_dense), 'k--', 'LineWidth', 2), hold on,
plot(x_dense,spline_fit(x_dense), 'LineWidth', 2)
scatter(x,f(x))
maxError = max( abs( f(x_dense)-spline_fit(x_dense) ) );
title(sprintf('chebyshev grid, N=%d, K=%d, max error = %.2g',N,K,maxError))


%%%%%%%%%%
% Now let's confirm we're doing the same thing as matlab.
x = [-1 -0.33 0.33 1]';
y = [0 1 2 0]';

spline_fit = InterpolatingSpline(x,y,4);
interpolant = griddedInterpolant(x,y,'spline');
figure
scatter(x,y), hold on,
plot(x_dense,spline_fit(x_dense), 'b', 'LineWidth', 4)
plot(x_dense,interpolant(x_dense), 'k--', 'LineWidth', 2)
