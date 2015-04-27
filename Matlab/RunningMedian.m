function [smoothVec] = RunningAverage( vec, smoothness )
%	Returns a running average of vec. Tries to handle end points well.

non_nan_indices = find( isnan(vec) == 0);

smoothVec = NaN(size(vec));
smoothVec(1) = vec(1);

smoothnessHalfLength = ceil((smoothness - 1)/2);

for index=2:length(vec)
	% Determine the largest width we're allowed for our running average=
	startDistance = index-1;
	endDistance = length(vec) - index;
	restrictionDistance = startDistance;
	if ( endDistance < restrictionDistance)
		restrictionDistance = endDistance;
	end
	
	if ( smoothnessHalfLength < restrictionDistance)
		restrictionDistance = smoothnessHalfLength;
	end
	
	indices = (index-restrictionDistance):(index+restrictionDistance);
	indices = intersect(indices, non_nan_indices);
	
	if ( length(indices) > 0 )
		smoothVec(index) = median(vec(indices));
	end
end



