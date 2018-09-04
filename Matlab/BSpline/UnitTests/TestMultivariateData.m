f = @(t) [sin(2*pi*t/10),cos(2*pi*t/10)];

t = (0:10)'; % observation points

% f_spline = BSpline(t,f(t));
f_spline = TensionSpline(t,f(t),0);

tq = linspace(min(t),max(t),10000)'; % evaluation points

figure
plot(tq,f_spline(tq))