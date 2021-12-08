% % AMS figure widths, given in picas, converted to points (1 pica=12 points)
% scaleFactor = 2;
% LoadFigureDefaults
% 
% FramesFolder = './DrifterSmoothingFrames';
% if exist(FramesFolder,'dir') == 0
% 	mkdir(FramesFolder);
% end
% 
% 
% shouldSaveFigures = 1;
% 
% % Drifter to highlight in the final plots
% choiceDrifter = 6;
% 
% drifters = open('../Matlab/sample_data/raw_rho1_drifters.mat');
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % We need a consistent lon0 for all drifters
% lon = drifters.lon{1};
% lon0 = min(lon)+(max(lon)-min(lon))/2;
% 
% drifter6 = GPSSmoothingSpline((drifters.date{6}-drifters.lastDeployment)*24*60*60,drifters.lat{6},drifters.lon{6},'shouldUseRobustFit',1,'lon0',lon0);
% % drifter7 = GPSSmoothingSpline((drifters.date{7}-drifters.lastDeployment)*24*60*60,drifters.lat{7},drifters.lon{7},'shouldUseRobustFit',1,'lon0',lon0);
% 
% pctmax = 1-1/100;
% zmin = 0;
% zmax = drifter6.noiseDistanceDistribution.locationOfCDFPercentile(pctmax);
% expectedVariance = drifter6.noiseDistanceDistribution.varianceInRange(zmin,zmax);
% 
% % grab the minimum value
% rmseRangedMin = sqrt(drifter6.expectedMeanSquareErrorInDistanceRange(zmin,zmax,expectedVariance)*drifter6.distribution.variance/2);
% rmseRangedMinLambda = drifter6.lambda;
% 
% drifter6.minimize( @(aSmoothingSpline) aSmoothingSpline.expectedMeanSquareErrorFromCV );
% rmseCVMin = sqrt(drifter6.expectedMeanSquareErrorFromCV/2);
% rmseCVMinLambda = drifter6.lambda;
% 
% 
% 
% nLambda = 200;
% lambdaRange = 10.^linspace(log10(rmseRangedMinLambda)-4,log10(rmseRangedMinLambda)+2,nLambda);
% 
% for iLambda = 1:nLambda
%     drifter6.lambda = lambdaRange(iLambda);
%     
%     mseCV(iLambda) = sqrt(drifter6.expectedMeanSquareErrorFromCV/2);
%     mseRanged(iLambda) = sqrt(drifter6.expectedMeanSquareErrorInDistanceRange(zmin,zmax,expectedVariance)*drifter6.distribution.variance/2);
% end
% 
% 
% t=linspace(min(drifter6.t),max(drifter6.t),length(drifter6.t)*10).';
% 
% s = 1/1000; % scale for x and y axis
% x0 = min(drifter6.x);
% y0 = drifter6.y(1);
% 
% datapointColor = 0*[1.0 1.0 1.0];
% datapointSize = 3.0*scaleFactor;
% outlierFaceColor = 1.0*[1.0 1.0 1.0];
% outlierEdgeColor = 0.2*[1.0 1.0 1.0];
% outlierSize = 5.5*scaleFactor;

figure('Position',[200 200 1280 720]);
set(gcf,'PaperPositionMode','auto');
set(gcf, 'Color', 'w');

iPlot = 0;
for iLambda = 1:nLambda
    iPlot = iPlot + 1;
    clf
    
    drifter6.lambda = lambdaRange(iLambda);
    
    subplot(1,2,1)
%     [xq,yq] = drifter7.xyAtTime(t);
%     p0 = plot(s*(xq-x0),s*(yq-y0), 'LineWidth', 1.0*scaleFactor, 'Color',0.8*[0.2 0.2 1.0]); hold on
    p1 = scatter(s*(drifter6.x(drifter6.outlierIndices)-x0),s*(drifter6.y(drifter6.outlierIndices)-y0),(outlierSize)^2, 'MarkerEdgeColor', outlierEdgeColor, 'MarkerFaceColor', outlierFaceColor); hold on
    p2 = scatter(s*(drifter6.x-x0),s*(drifter6.y-y0),(datapointSize)^2,'filled', 'MarkerEdgeColor', datapointColor, 'MarkerFaceColor', datapointColor);
    
    [xq,yq] = drifter6.xyAtTime(t);
    
    p3 = plot(s*(xq-x0),s*(yq-y0), 'LineWidth', 1.0*scaleFactor, 'Color',0.4*[1.0 1.0 1.0]);
    
    xlabel('x (km)', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
    ylabel('y (km)', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
    axis equal
    xlim([-0.5 3.5])
    ylim([-0.5 5.0])
    
     leg = legend([p3 p2 p1],'spline fit', 'data', 'outlier');
    leg.FontSize = figure_legend_size;
    leg.FontName = figure_font;
    leg.Location = 'northwest';
    subplot(1,2,2)
    
    plot(lambdaRange,mseCV,'LineWidth', 2*scaleFactor), xlog, ylog
    hold on
    plot(lambdaRange,mseRanged,'LineWidth', 2*scaleFactor)
    scatter( rmseRangedMinLambda, rmseRangedMin, (10*scaleFactor)^2, 'filled')
    scatter( rmseCVMinLambda, rmseCVMin, (10*scaleFactor)^2, 'filled')
    plot([lambdaRange(iLambda) lambdaRange(iLambda)],[100 500],'k','LineWidth', 2*scaleFactor)
    ylim([100 500])
    xlim([min(lambdaRange) max(lambdaRange)])
    ax = gca;
    ax.YAxisLocation = 'right';
    xlabel('\lambda', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
    leg = legend('CV MSE','Ranged MSE');
    leg.FontSize = figure_legend_size;
    leg.FontName = figure_font;

    packfig(1,2)
%     tightfig
    
    print('-depsc2', sprintf('%s/t_%03d', FramesFolder,iPlot) );
    if lambdaRange(iLambda) <= rmseRangedMinLambda && lambdaRange(iLambda+1) > rmseRangedMinLambda
        for i=1:47
            iPlot = iPlot + 1;
            print('-depsc2', sprintf('%s/t_%03d', FramesFolder,iPlot) );
        end
    elseif lambdaRange(iLambda) <= rmseCVMinLambda && lambdaRange(iLambda+1) > rmseCVMinLambda
        for i=1:47
            iPlot = iPlot + 1;
            print('-depsc2', sprintf('%s/t_%03d', FramesFolder,iPlot) );
        end
    end
end
