% AMS figure widths, given in picas, converted to points (1 pica=12 points)
scaleFactor = 2;
LoadFigureDefaults
addpath('support')

shouldSaveFigures = 1;

% Drifter to highlight in the final plots
choiceDrifter = 6;

drifters = open('../Matlab/sample_data/raw_rho1_drifters.mat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We need a consistent lon0 for all drifters
lon = drifters.lon{1};
lon0 = min(lon)+(max(lon)-min(lon))/2;

drifter6 = GPSSmoothingSpline((drifters.date{6}-drifters.lastDeployment)*24*60*60,drifters.lat{6},drifters.lon{6},'shouldUseRobustFit',1,'lon0',lon0);
drifter7 = GPSSmoothingSpline((drifters.date{7}-drifters.lastDeployment)*24*60*60,drifters.lat{7},drifters.lon{7},'shouldUseRobustFit',1,'lon0',lon0);


s = 1/1000; % scale for x and y axis
x0 = min(drifter6.x);
y0 = drifter7.y(1);

datapointColor = 0*[1.0 1.0 1.0];
datapointSize = 3.0;

figure
subplot(1,2,1)
plot(s*(drifter6.x-x0),s*(drifter6.y-y0)), hold on
scatter(s*(drifter6.x-x0),s*(drifter6.y-y0),(datapointSize)^2,'filled', 'MarkerEdgeColor', datapointColor, 'MarkerFaceColor', datapointColor)


title('drifter 6', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
xlabel('x (km)', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
ylabel('y (km)', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
axis equal
xlim([-0.5 3.5])
ylim([-0.5 4.5])

subplot(1,2,2)
plot(s*(drifter7.x-x0),s*(drifter7.y-y0)), hold on
scatter(s*(drifter7.x-x0),s*(drifter7.y-y0),(datapointSize)^2,'filled', 'MarkerEdgeColor', datapointColor, 'MarkerFaceColor', datapointColor)


title('drifter 7', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
xlabel('x (km)', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
% ylabel('y (km)', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
axis equal
xlim([-0.5 3.5])
ylim([-0.5 4.5])

packfig(1,2)
tightfig

print('-depsc2', 'figures/GPSRawData.eps')