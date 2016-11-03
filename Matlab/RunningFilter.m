% Computes a running average across ensembles and reports standard error
% size(vec) = [n m] where n is time, and m is the ensemble

function [y,yerr] = RunningFilter( x, smoothness, filter )
%	Returns a running average of vec. Tries to handle end points well.

n = size(x,1);
m = size(x,2);

y = NaN(n,1);
yerr = NaN(n,1);

smoothnessHalfLength = ceil((smoothness - 1)/2);

for index=1:n
    restrictionDistance = min([smoothnessHalfLength; index-1; n-index]);
    indices = (index-restrictionDistance):(index+restrictionDistance);
    
    if ( ~isempty(indices) )
        vals = reshape(x(indices,:),[length(indices)*m 1]);
        if strcmp(filter,'mean')
            y(index) = mean(vals);
        elseif strcmp(filter,'median')
            y(index) = median(vals);
        end
        yerr(index) = std(vals)/sqrt(length(indices)*m);
    end
end



