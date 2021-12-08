FramesFolder = './SmoothingParameterTrueFunctionFrames';
if exist(FramesFolder,'dir') == 0
	mkdir(FramesFolder);
end

L = 10;
f = @(x) sin(2*pi*x/L);

N = 30;
x = linspace(0,L,N)';

sigma = 0.2;
y = f(x) + sigma*randn(N,1);

x_dense = linspace(0,L,10*N)';

spline = SmoothingSpline(x,y,NormalDistribution(sigma));
lambda_optimal = spline.lambda;

lambdaRange = 10.^linspace(-6,6,96);

figure('Position',[200 200 1280/2 720/2]);
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Color', 'w');

for i=1:length(lambdaRange)
    clf
    spline.lambda = lambdaRange(i);
    plot(x_dense,f(x_dense)), hold on
    plot(x_dense,spline(x_dense),'LineWidth',2,'Color',[0.9290    0.6940    0.1250]), hold on
    
    scatter(x,y,10^2,'filled','k')
    ylim([-1.35 1.35])
    text( 7, 1, sprintf('\\lambda = %.1e',spline.lambda),'Interpreter','tex', 'LineWidth', 2, 'FontName', 'Helvetica', 'FontSize', 24, 'BackgroundColor', 'white')
    
    print('-depsc2', sprintf('%s/t_%03d', FramesFolder,i) );
    print('-depsc2', sprintf('%s/t_%03d', FramesFolder,2*length(lambdaRange)-i+1) );
end

clf
spline.lambda = lambda_optimal;
plot(x_dense,f(x_dense)), hold on
plot(x_dense,spline(x_dense),'LineWidth',2,'Color',[0.9290    0.6940    0.1250]), hold on
    
scatter(x,y,10^2,'filled','k')
ylim([-1.35 1.35])
text( 7, 1, sprintf('\\lambda = %.1e',spline.lambda),'Interpreter','tex', 'LineWidth', 2, 'FontName', 'Helvetica', 'FontSize', 24, 'BackgroundColor', 'white')

print('-depsc2', sprintf('%s/expected_least_square_error', FramesFolder) );

% legend('true function', 'noisy data', 'tension spline fit')

% spline.lambda = 0;