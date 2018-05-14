function [m_x,Cm_x,B] = bspline_fit_no_tension(t,x,dx,S,t_knot,weight_function)
% bspline_fit_no_tension    Find the maximum likelihood fit
%
% t         independent variable (time), length N
% x         observations at time t_i, length N
% dx        error of observation x, length N
% S         degree of spline (e.g., 3 denotes cubic), scalar
% t_knot    knot points for the splines
% weight_function   used for iteratively reweighted least-squares
%
% 
% B is a matrix containing the M B-splines, at the N locations t_i, for
% K=S+1 derivatives. The matrix is NxMxK.
%
% m_x       The coefficients for the model fit, length M.
% Cm_x      The error in the coefficients, size MxM

if (length(t) ~= length(x) )
    disp('The time series are not consistent lengths');
    return;
end

B = bspline(t,t_knot,S+1);
X = squeeze(B(:,:,1));
% N = size(X,1); % also size(X,1);
% M = size(X,2); % number of splines

Wx = diag(1./(dx.^2));

dbstop if warning
[m_x,Cm_x] = ComputeSolution( X, Wx, x );

error_x_previous = dx;
rel_error = 1.0;
repeats = 1;
while (rel_error > 0.01)
    dx2 = weight_function(X*m_x - x);
    
    Wx = diag(1./(dx2));
    
    [m_x,Cm_x] = ComputeSolution( X, Wx, x );
    
    rel_error = max( (dx2-error_x_previous)./dx2 );
    error_x_previous=dx2;
    repeats = repeats+1;
    
    if (repeats == 100)
        disp('Failed to converge after 100 iterations.');
        break;
    end
end

end


function [m_x, Cm_x] = ComputeSolution( X, Wx, x )
% X matrix:
% Rows are the N observations
% Columns are the M splines
% Wx is NxN
% x is Nx1

% set up inverse matrices
E_x = X'*Wx*X; % MxM

m_x = E_x\(X'*Wx*x);

Cm_x = inv(E_x);
end

