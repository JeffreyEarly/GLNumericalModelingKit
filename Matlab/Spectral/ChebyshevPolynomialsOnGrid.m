function [varargout] = ChebyshevPolynomialsOnGrid( x, N_polys )
%% Chebyshev Polynomials on Grid
% Compute the the first N Chebyshev polynomials and their derivatives for
% an arbitrary grid x.
% 
% x_norm = ChebyshevPolynomialsOnGrid( x ) with exactly one argument, x,
% returns the x normalized to its typical [-1,1] values.
%
% T = ChebyshevPolynomialsOnGrid( x, N_polys ) returns the first N_poly
% Chebyshev polynomials for an arbitrary grid x.
%
% [T, T_x, T_xx,...] = ChebyshevPolynomialsOnGrid( x, N_polys ) returns the
% first N_poly Chebyshev polynomials and their derivatives for an arbitrary
% grid x.
%
% The returned matrices T, T_xx, etc are size(T) = [length(x) N_polys],
% i.e., the polynomials are given column-wise.

    N_points = length(x);
    
    % These are the normalized coordinates for Chebyshev polynomials.
    L = max(x)-min(x);
    x_norm = (2/L)*(x-min(x)) - 1;
    t = acos(x_norm);
    
    % if there's only one input argument, we just return x_norm
    if nargin == 1
        varargout{1} = x_norm;
        return;
    else
        if N_polys < 4
           disp('You must request at least four polynomials. Fixing that for you.')
           N_polys = 4;
        end
    end
    
    N_diff = nargout-1;
    varargout = cell(1,nargout);
    
    % It's easy to create the polynomials, they're stretched cosines!
    T = zeros(N_points,N_polys);
    for iPoly=0:(N_polys-1)
       T(:,iPoly+1) = cos(iPoly*t);
    end
    
    varargout{1} = T;
    
    % Now use the recursion formula to compute derivates of polynomials.
    for n=1:N_diff
        T = varargout{n};
        T_x = zeros(size(T));
        T_x(:,2) = T(:,1);
        T_x(:,3) = 2*2*T(:,2);
        for j=4:N_polys
            m = j-1;
            T_x(:,j) = (m/(m-2))*T_x(:,j-2) + 2*m*T(:,j-1);
        end
        varargout{n+1} = (2/L)*T_x;
    end
    
end