% The eigenvalues and modes are ordered from highest variances, to lowest.
% Each column of the mode matrix contains an eigenmode, e.g., mode(:,1) contains the mode
% the most variance.
%
% The score matrix contains the weights of each mode necessary to reconstruct the observations.
% So: recoveredObservations = (score(:,1:3) * modes(:,1:3)' ) + repmat(obsMean,[n 1]);
function [modes, eigenvalues, score, obsMean] = PrincipalComponentsOfObservations( observations )

% Compute the mean and standard deviation.
obsMean = mean(observations);
obsStd = std(observations);

% n = number of observations
% m = number of depth variables
[n m] = size(observations);

% Subtract the mean, and divide by the standard deviation to normalize the data
normalizedObservations = (observations - repmat(obsMean,[n 1])) ./ repmat(obsStd,[n 1]);

% Find the eigenfunctions of the covariance matrix.
% The rows of the modes matrix represent the individual modes
% Matrices are m x m
[modes eigenvalueMatrix] = eig(cov(normalizedObservations));

[eigenvalues, permutation] = sort(real(diag(eigenvalueMatrix)),1,'descend');
modes=modes(:,permutation);

% This adds the variance back in to each mode.
modes = modes .* repmat(sqrt(obsStd'),[1 m]);

% Score is now #obs by #depth points, (nxm).(mxm) = nxm
% Each row represents an observation, each column represents the weighting of the particular mode.
score = normalizedObservations * modes;