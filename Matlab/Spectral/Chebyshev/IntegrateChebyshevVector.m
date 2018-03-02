function v_p = IntegrateChebyshevVector(v)
% Taken from cumsum as part of chebfun
%
% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% integration target
n = length(v);
v_p = zeros(n+1,1);

% zero-pad coefficients
v = reshape(v,[],1);
v = [v; zeros(2,1)];

v_p(3:n+1) = (v(2:n) - v(4:n+2)) ./ (2*(2:n).'); % Compute b_(2) ... b_(n+1):
v_p(2) = v(1) - v(3)/2;        % Compute b_1

t = ones(1, n);
t(2:2:end) = -1;
v_p(1) = t*v_p(2:end);             % Compute b_0 (satisfies f(-1) = 0)
end