D = importdata('motionless_garmin_edge_705.txt','\t');

lat = D.data(:,1);
lon = D.data(:,2);
t = datetime(D.textdata(6:end,2));

indices = 1:5:5*1200;
spline = GPSTensionSpline(t(indices),lat(indices),lon(indices));

tq = linspace(min(spline.t),max(spline.t),10*length(spline.t))';

[x,y] = spline.xyAtTime(tq);
figure
scatter(spline.t,spline.x), hold on
plot(tq,x)