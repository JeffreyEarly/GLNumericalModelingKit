K = 3; % order of spline
D = K-1; % number of derivates to return
t = (0:10)'; % observation points

% increase the multiplicity of the end knots for higher order splines
t_knot = [repmat(t(1),K-1,1); t; repmat(t(end),K-1,1)];

nSplines = length(t)+K-2;

% coefficients for the bsplines---set all of them to zero for now.
m = zeros(nSplines,1);

% initialize the BSpline class
B = BSpline(K,t_knot,m);

% make a really fine grid to draw the splines on
tq = linspace(min(t),max(t),1000)';

figure
subplot(2,1,1)
for iSpline = 1:nSplines
   m = zeros(nSplines,1);
   m(iSpline) = 1;
   B.m = m;
   plot(tq,B(tq),'LineWidth',2), hold on
end
title('b-splines')
subplot(2,1,2)
for iSpline = 1:nSplines
   m = zeros(nSplines,1);
   m(iSpline) = 1;
   B.m = m;
   plot(tq,B(tq,1),'LineWidth',2), hold on
end
title('1st derivative')
print('-depsc2', '../figures/bspline.eps')