classdef ConstrainedSpline < BSpline
    %ConstrainedSpline Summary of this class goes here
    %   2 argument initialization
    %       f = ConstrainedSpline(t,x,K,t_knot,constraints);
    %   where
    %       t               array of values for the independent axis
    %       x               array of values for the dependent axis 
    %       K               spline order
    %       t_knot          array of knot points
    %       constraints     constraints = struct('t',[],'D',[]);
    %       f               cubic spline interpolant
    % 
    % constraints are such that f^(D)(t)=0
    
    
    properties (Access = public)
        distribution
        t
        x
        
        CmInv
        X
        W
    end
    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Initialization
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function self = ConstrainedSpline(t,x,K,t_knot,distribution,constraints)
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
            
            if isempty(distribution)
                distribution = NormalDistribution(1);
            end
            
            nl = find(t_knot <= t_knot(1),1,'last');
            nr = length(t_knot)- find(t_knot == t_knot(end),1,'first')+1;
            t_knot = [repmat(t_knot(1),K-nl,1); t_knot; repmat(t_knot(end),K-nr,1)];
            
            % Find the spline coefficients
            X = BSpline.Spline( t, t_knot, K, 0 );
            
            % Deal with constraints, if any
            if isempty(constraints)
                F=[];
            else
                M = size(X,2); % number of splines
                NC = length(constraints.t); % number of constraints
                if length(constraints.D) ~= NC
                    error('t and D must have the same length in the constraints structure.');
                end
                F = zeros(NC,M);
                Xc = BSpline.Spline( constraints.t, t_knot, K, K-1 );
                for i=1:NC
                    F(i,:) = squeeze(Xc(i,:,constraints.D(i)+1));
                end

            end
            if isa(distribution,'NormalDistribution')
                [m,CmInv] = ConstrainedSpline.ConstrainedSolution(X,x,F,distribution.sigma);
                W = 1/(distribution.sigma0)^2;
            else
                [m,CmInv,W] = ConstrainedSpline.IteratedLeastSquaresConstrainedSolution(X,x,F,sqrt(distribution.variance),[],[],distribution.w);
            end
            
            self@BSpline(K,t_knot,m);
            self.distribution = distribution;
            self.t = t;
            self.x = x;
            self.CmInv = CmInv;
            self.X = X;
            self.W = W;
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Smoothing matrix and covariance matrix
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function S = smoothingMatrix(self)
            % The smoothing matrix S takes the observations and maps them
            % onto the estimated true values.
            S = (self.X*(self.CmInv\(self.X.')))*self.W;
        end
    end
    
    
    methods (Static)
        
        function tc = MinimumConstraintPoints(t_knot,K,T)
            %% MinimumConstraintPoints
            % Given some set of knot points for a spline of order K, if you
            % want to universally constrain that spline at degree T, this
            % function returns a set of constraint locations (in t) that
            % let you do that, assuming repeat knot points at the beginning
            % and end.
            
            t = unique(t_knot);
            D = K-1-T; % 0 if we're constraining at the same order
            if mod(D,2) == 0
                ts = t(1) + (t(2)-t(1))/(D/2 + 2)*(1:(D/2+1)).';
                te = t(end-1) + (t(end)-t(end-1))/(D/2 + 2)*(1:(D/2+1)).';
                ti = t(2:end-2) + diff(t(2:end-1))/2;
                tc = cat(1,ts,ti,te);
            else
                ts = t(1) + (t(2)-t(1))/((D-1)/2+1)*(0:((D-1)/2)).';
                te = t(end) - (t(end)-t(end-1))/((D-1)/2 + 1)*(0:((D-1)/2)).';
                ti = t(2:end-1);
                tc = cat(1,ts,ti,te);
            end
        end
        
        function [m,CmInv] = ConstrainedSolution(X,x,F,sigma,XWX,XWx)
            %% ConstrainedSolution
            %
            % Returns the spline fit solution with constraints, given by
            % the F matrix. It also supports IRLS with non-constant sigma.
            %
            %
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
            
            if ~exist('sigma','var') || isempty(sigma)
                sigma = 1;
            end
            
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
            F_x = XWx;
            
            if NC > 0
                E_x = cat(1,E_x,F); % (M+NC)xM
                E_x = cat(2,E_x,cat(1,F',zeros(NC)));
                F_x = cat(1,F_x,zeros(NC,1));
                m_x = E_x\F_x;
                m = m_x(1:size(X,2));
            else
                m = E_x\F_x;
            end
            
            
            CmInv = XWX;
        end
        
        function [m,CmInv,W] = IteratedLeastSquaresConstrainedSolution(X,x,F,sigma,XWX,XWx,w,t,rho)
            % Same calling sequence as the TensionSolution function, but
            % also includes the weight factor, w
            if exist('t','var') && exist('rho','var') && ~isempty(rho)
                rho_t = rho(t - t.');
                if length(sigma) == 1
                    sigma = ones(size(x))*sigma;
                end
                Sigma2 = (sigma * sigma.') .* rho_t;
                W = inv(Sigma2);
                m = ConstrainedSpline.ConstrainedSolution(X,x,F,W,XWX,XWx);
            else
                m = ConstrainedSpline.ConstrainedSolution(X,x,F,sigma,XWX,XWx);
            end
            
            
            sigma2_previous = sigma.*sigma;
            rel_error = 1.0;
            repeats = 1;
            while (rel_error > 0.01)
                sigma_w2 = w(X*m - x);
                
                if exist('rho_t','var')
                    Sigma2 = (sqrt(sigma_w2) * sqrt(sigma_w2).') .* rho_t;
                    W = inv(Sigma2);
                else
                    W = diag(1./sigma_w2);
                end
                
                [m,CmInv] = ConstrainedSpline.ConstrainedSolution(X,x,F,W,XWX,XWx);
                
                rel_error = max( abs((sigma_w2-sigma2_previous)./sigma_w2) );
                sigma2_previous=sigma_w2;
                repeats = repeats+1;
                
                if (repeats == 150)
                    disp('Failed to converge after 150 iterations.');
                    break;
                end
            end
            
        end
  
    end
end


