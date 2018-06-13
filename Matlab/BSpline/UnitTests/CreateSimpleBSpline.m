K = 3; % order of spline
D = K-1; % number of derivates to return
t = (0:10)'; % observation points
t_knot = BSpline.KnotPointsForPoints(t,K);

tq = linspace(min(t),max(t),10000)';

B = BSpline.Spline( tq, t_knot, K, D );

figure
for iPlot = 1:size(B,3)
   subplot(size(B,3),1,iPlot)
   plot(tq,squeeze(B(:,:,iPlot)))
end
