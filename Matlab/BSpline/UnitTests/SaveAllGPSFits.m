drifters = open('raw_rho2_drifters.mat');
nDrifters = length(drifters.date);

% figure
% for iDrifter=1:length(drifters.date)
%     t_drifter = (drifters.date{iDrifter}-drifters.lastDeployment)*24*60*60;
%     scatter(t_drifter/3600,1*ones(size(t_drifter))), hold on;
% end


offsets = 0:14;
distances = zeros(length(offsets),length(drifters.date));
for iOffset = 1:length(offsets)
    offset = offsets(iOffset);
    for iDrifter=1:length(drifters.date)
        t_zero = drifters.lastDeployment*24*60 - offset;
        t_drifter = drifters.date{iDrifter}*24*60 - t_zero;
        t_drifter_mod = mod(t_drifter,30);
        t_drifter_mod(t_drifter_mod>15)=t_drifter_mod(t_drifter_mod>15)-15;
        distances(iOffset,iDrifter) = median(t_drifter_mod);
    end
end

% offset of 2 appears to be minimum distances.
[val,idx] = min(median(distances,2));

t_zero = drifters.lastDeployment*24*60*60 - offsets(idx)*60;
t_end = drifters.firstRetrieval*24*60*60 - t_zero;
t_interp = (0:1800:t_end).';

x_drifters = nan(length(t_interp),nDrifters);
y_drifters = nan(length(t_interp),nDrifters);
n_eff_x = nan(nDrifters,1);
n_eff_y = nan(nDrifters,1);

for iDrifter=1:nDrifters
    fprintf('..%d',iDrifter);
    t_drifter = drifters.date{iDrifter}*24*60*60 - t_zero;
    nonExtrapIndices = t_interp >= t_drifter(1);
    if sum(nonExtrapIndices) < length(t_interp)
        fprintf('only using %d of %d indices to avoid extrapolation\n',sum(nonExtrapIndices),length(t_interp));
    end
    
    if iDrifter == 1
        spline = GPSSmoothingSpline(t_drifter,drifters.lat{iDrifter},drifters.lon{iDrifter},'shouldUseRobustFit',1);
        [x,y] = spline.xyAtTime(t_interp);
        
        lon0 = spline.lon0;
        x0 = min(x);
        y0 = min(y);
        x = x-x0;
        y = y-y0;
    else
        spline = GPSSmoothingSpline(t_drifter,drifters.lat{iDrifter},drifters.lon{iDrifter},'shouldUseRobustFit',1,'lon0',lon0,'x0',x0,'y0',y0);
        [x,y] = spline.xyAtTime(t_interp(nonExtrapIndices));
    end
    
    [n_x,n_y] = spline.effectiveSampleSize;
    x_drifters(nonExtrapIndices,iDrifter) = x;
    y_drifters(nonExtrapIndices,iDrifter) = y;
    n_eff_x(iDrifter) = n_x;
    n_eff_y(iDrifter) = n_y;
end
fprintf('\n');

x = x_drifters;
y = y_drifters;
t = t_interp;
lat0 = TransverseMercatorToLatitudeLongitude(x0,y0,lon0,0.9996);

save('smoothedGriddedRho2Drifters.mat','t','x','y','x0','y0','lon0','lat0','n_eff_x','n_eff_y');

figure
plot(x,y), axis equal

% figure
% scatter(spline.x-x0,spline.y-y0,(2.5)^2,'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k'), hold on
% scatter(x,y,(2.5)^2,'filled', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r')