f = @(x) x.^2 + x - 1;
f = @(x) sin(4*pi*x);

xMin = -1;
xMax = 1;

[xLobatto, f_cheb] = ProjectOntoChebyshevPolynomialsWithTolerance(f, [xMin xMax], 1e-16);

fx_cheb = DifferentiateChebyshevVector(f_cheb);

roots = FindRootsFromChebyshevVector(f_cheb, xLobatto);
extrema = FindRootsFromChebyshevVector(fx_cheb(1:end-1), xLobatto);

x = linspace(xMin,xMax,1000)';
figure, plot(x,f(x)), hold on
scatter(roots,f(roots))
scatter(extrema,f(extrema))
