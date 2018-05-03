function [AC, DOF] = Autocorrelation(x,maxlag)
% Autocorrelation
% size(v) = [N M] where N is length in time, M are reps.
%
% Returns the unbiased estimator.
%
% We do *not* remove the mean from x. The expected mean may be zero, or
% some other value.

% v = v-mean(v);
% sigma2 = std(v).^2;
sigma2 = mean(x.*x);
AC = zeros(size(x));
DOF = length(x) - (0:maxlag)';
for lag=0:maxlag
    x_shift = circshift(x,-lag,1);
    v2 = x.*x_shift;
     AC(lag+1,:) = mean(v2(1:(length(x)-lag))); % unbiased estimator
%    AC(lag+1) = sum(v2(1:(length(v)-lag)))/DOF(lag+1); % Proper accounting
end
AC = AC(1:(maxlag+1),:)./sigma2;