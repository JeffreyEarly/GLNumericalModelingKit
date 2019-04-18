drifters = open('raw_rho1_drifters.mat');

iDrifter = 6;

t_drifter = (drifters.date{iDrifter}-drifters.lastDeployment)*24*60*60;
spline = GPSSmoothingSpline(t_drifter,drifters.lat{iDrifter},drifters.lon{iDrifter},'shouldUseRobustFit',1);
tq = linspace(min(t_drifter),max(t_drifter),10*length(t_drifter))';
[x,y] = spline.xyAtTime(tq);
figure
scatter(spline.x,spline.y,(2.5)^2,'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k'), hold on
plot(x,y),axis equal
