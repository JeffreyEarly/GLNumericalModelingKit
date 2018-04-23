function [AC, DOF] = Autocorrelation(v,maxlag)
% Autocorrelation
% size(v) = [N M] where N is length in time, M are reps

% v = v-mean(v);
% sigma2 = std(v).^2;
sigma2 = mean(v.*v);
AC = zeros(size(v));
DOF = length(v) - (0:maxlag)';
for lag=0:maxlag
    v_shift = circshift(v,-lag,1);
    v2 = v.*v_shift;
     AC(lag+1,:) = mean(v2(1:(length(v)-lag))); % Crude estimate
%    AC(lag+1) = sum(v2(1:(length(v)-lag)))/DOF(lag+1); % Proper accounting
end
AC = AC(1:(maxlag+1),:)./sigma2;