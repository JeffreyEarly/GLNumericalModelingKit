scaleFactor = 2;
LoadFigureDefaults

sigma = 8.5;
nu = 4.5;
alpha = 0.5;

dist1 = StudentTDistribution(sigma,nu);
dist2 = StudentTDistribution(10*sigma,nu);
totalDist = AddedDistribution(alpha,dist2,dist1);

x = linspace(-150,150,1000)';

FigureSize = [50 50 figure_width_1col+4 150*scaleFactor];
fig1 = figure('Units', 'points', 'Position', FigureSize);
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Color', 'w');
fig1.PaperUnits = 'points';
fig1.PaperPosition = FigureSize;
fig1.PaperSize = [FigureSize(3) FigureSize(4)];

plot(x,dist1.pdf(x),'LineWidth',2*scaleFactor), hold on
plot(x,dist2.pdf(x),'LineWidth',2*scaleFactor)
plot(x,totalDist.pdf(x),'k','LineWidth',2*scaleFactor)

legend('GPS','Outliers','Total')

xlabel('meters', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
set( gca, 'FontSize', figure_axis_tick_size);

print('-depsc2', 'robust_distribution.eps')