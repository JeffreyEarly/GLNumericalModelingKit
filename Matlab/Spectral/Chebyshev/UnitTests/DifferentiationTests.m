f = @(x) x.^2 + x - 2;
f_x = @(x) 2*x + 1;


xMin = -1;
xMax = 1;
L = xMax-xMin;
n = 10;

xLobatto = (L/2)*( cos(((0:n-1)')*pi/(n-1)) + 1) + xMin;

f_cheb = fct(f(xLobatto));
fx_cheb = DifferentiateChebyshevVector(f_cheb);

F_x = ifct(fx_cheb);
Dx = ChebyshevInterpolationDerivative(n);

figure
subplot(3,1,1)
plot(xLobatto,F_x), hold on
scatter(xLobatto,f_x(xLobatto))
title('fct,diff, ifct')
subplot(3,1,2)
plot(xLobatto,Dx*f(xLobatto)), hold on
scatter(xLobatto,Dx*f(xLobatto))
title('diff')

fx_int_cheb = IntegrateChebyshevVector(fx_cheb);
F = ifct(fx_int_cheb(1:n));

subplot(3,1,3)
plot(xLobatto,F), hold on
scatter(xLobatto,f(xLobatto)+2)
title('int')


