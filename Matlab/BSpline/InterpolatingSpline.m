classdef InterpolatingSpline < BSpline
    %InterpolatingSpline Summary of this class goes here
    %   2 argument initialization
    %       f = InterpolatingSpline(t,x);
    %   where
    %       t       array of values for the independent axis
    %       x       array of values for the dependent axis 
    %       f       cubic spline interpolant
    % 
    %   3 argument initialization
    %       f = InterpolatingSpline(t,x,K);
    %   where
    %       K       order of the spline
    
    
    properties (Access = public)

    end
    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Initialization
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function self = InterpolatingSpline(t,x,varargin)
            % Make sure (t,x) are the right shapes, and see how many
            % dimensions of x we are given.
            t = reshape(t,[],1);
            N = length(t);
            x = reshape(x,[],1);
            if length(x) ~= N
                error('x and t must have the same length.');
            end
            
            if any(isnan(x)) || any(isnan(t))
                error('Some of the data contain NaNs!')
            end
            
            % Parse any extra input options.
            nargs = length(varargin);
            if nargs == 1
                K = varargin{1};
            elseif nargs == 0
                K = 4;
            elseif mod(nargs,2) == 0
                for k = 1:2:length(varargin)
                    if strcmp(varargin{k}, 'K')
                        K = varargin{k+1};
                    elseif strcmp(varargin{k}, 'S')
                        K = varargin{k+1}+1;
                    end
                end
            else
                error('Arguments must be given as name/value pairs');
            end
            
            % create knot points for a canonical interpolating spline
            t_knot = InterpolatingSpline.KnotPointsForPoints(t,K);
            
            % Find the spline coefficients
            X = BSpline.Spline( t, t_knot, K, 0 );
            m = X\x;
            
            self@BSpline(K,t_knot,m);
            
        end

    end
    
    
    methods (Static)
        
        function [m,CmInv] = ConstrainedSolution(X,sigma,x,F,XWX,XWx)
            % N     # of observations
            % M     # of splines
            % NC    # of constraints
            %
            % inputs:
            % X         splines on the observation grid, NxM
            % sigma     errors of observations, either a scalar, Nx1, OR if
            %           size(sigma)=[N N], then we assume it's the weight
            %           matrix W
            % x         observations (Nx1)
            % XWX       (optional) precomputed matrix X'*Wx*X
            % XWx       (optional) precomputed matrix X'*Wx*x
            %
            % output:
            % m         coefficients of the splines, Mx1
            % CmInv     Inverse of the covariance of coefficients, MxM
            %
            % X is NxM
            % W is NxN
            % x is Nx1
            % F is NCxM
            N = length(x);
            NC = size(F,1);
            
            if ~exist('XWX','var') || isempty(XWX)
                if size(sigma,1) == N && size(sigma,2) == N
                    XWX = X'*sigma*X;
                elseif length(sigma) == 1
                    XWX = X'*X/(sigma*sigma);
                elseif length(sigma) == N
                    XWX = X'*diag(1./(sigma.^2))*X; % MxM
                else
                    error('sigma must have the same length as x and t.');
                end
            end
            
            if ~exist('XWx','var') || isempty(XWx)
                if size(sigma,1) == N && size(sigma,2) == N
                    XWx = X'*sigma*x;
                elseif length(sigma) == 1
                    XWx = X'*x/(sigma*sigma);
                elseif length(sigma) == N
                    XWx = X'*diag(1./(sigma.^2))*x; % (MxN * NxN * Nx1) = Mx1
                else
                    error('sigma must have the same length as x and t.');
                end
            end
            
            % set up inverse matrices
            E_x = XWX; % MxM
            E_x = cat(1,E_x,F); % (M+NC)xM
            E_x = cat(2,E_x,cat(1,F',zeros(NC)));
            
            F_x = XWx;
            F_x = cat(1,F_x,zeros(NC,1));
            
            m_x = E_x\F_x;
            m = m_x(1:size(X,2));
            
            CmInv = XWX;
        end
        
        function t_knot = KnotPointsForPoints( t, K, DF )
            %% KnotPointsForPoints
            %
            % Returns the natural knot points for splines of order K for
            % observations at time points t. These knot points are place so
            % that each spline has exactly enough support for a fully
            % constrained system.
            %
            % It defaults to 1 degree of freedom, but you can override this
            % by setting DF to some nonnegative integer value.
            %
            % This matches the definitions in interpolation with tension
            % splines paper.
            
            if nargin < 3
                DF = 1;
            end
            
            if (DF < 1 || mod(DF,1) ~= 0)
                disp('DF must be a non-negative integer');
                return;
            end
            
            t = sort(t);
            
            t = [t(1); t(1+DF:DF:end-DF); t(end)];
            
            if mod(K,2) == 1
                % Odd spline order, so knots go in between points.
                dt = diff(t);
                
                % This gives us N+1 knot points
                t_knot = [t(1); t(1:end-1)+dt/2; t(end)];
                
                % Now remove start and end knots
                for i=1:((K-1)/2)
                    t_knot(2) = [];
                    t_knot(end-1) = [];
                end
                
            else
                t_knot = t;
                
                % Now remove start and end knots
                for i=1:((K-2)/2)
                    t_knot(2) = [];
                    t_knot(end-1) = [];
                end
                
            end
            
            % Now we increase the multiplicity of the knot points at the beginning and
            % the end of the interval so that the splines do not extend past the end
            % points.
            t_knot = [repmat(t_knot(1),K-1,1); t_knot; repmat(t_knot(end),K-1,1)];
        end
        
        function t_knot = KnotPointsForSplines( t, K, N_splines )
            %% KnotPointsForSplines
            %
            % Creates knot points for the data such that there will be
            % N_splines, each of which have approximately the same number
            % of observation points in support.
            %
            % Very crude, quick way of doing this.
            
            if N_splines < K
                N_splines = K;
            end
            
            t_knot = BSpline.KnotPointsForPoints(t, K, ceil(length(t)/N_splines));
        end
  
    end
end


