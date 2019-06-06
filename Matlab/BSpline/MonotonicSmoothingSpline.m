classdef MonotonicSmoothingSpline < SmoothingSpline
    %SmoothingSpline Fit noisy data with a tensioned interpolating spline
    %   3 argument initialization
    %       f = SmoothingSpline(t,x,sigma);
    %   where
    %       t               array of values for the independent axis
    %       x               array of values for the dependent axis 
    %       distribution    distribution of the noise
    %   returns
    %       f       spline interpolant
    %
    %   SmoothingSpline takes a number of optional input argument pairs.
    %
    %   The distribution must be a subclass of Distribution class.
    %
    %   'lambda' lambda is the tension parameter, and can be given directly
    %   as a numeric value, or can be a function that takes this
    %   SmoothingSpline object as an argument, and returns a numeric value.
    %
    % 
    
    properties
        x
        t
        distribution
        
        T           % degree at which tension is applied
        lambda      % tension parameter
        isConstrained % indicates whether or not lambda was so big that the solution is just a constrained solution
        mu          % mean value of the tension variable
        knot_dof    % knot dofs
        
        covariance  % computed from the given distribution, this is the covariance structure of the observations. It may be a scalar, vector, or matrix.
        
        variableCache % structure storing several cached variables, useful for quick tension spline computation.
        
        sigma       % initial weight (given as normal std dev.)
        
        
        % These have no consequence to the fit, and are only populated if
        % 'shouldEstimateOutlierDistribution' is set to 1.
        outlierDistribution = []
        alpha = []
        lambdaAtFullTension % what was the value of lambda at full tension
        sigmaAtFullTension % what was the set of 'sigmas' produced from the full tension solution
        
        % These have no consequence to the fit, but can be useful for
        % diagnostics.
        outlierIndices = [] 
        outlierThreshold % set to a distance with < 1/10000 odds
    end
    
    properties (Dependent)
        X % Splines at observation points, NxM
        W % Weight matrix at observation points, NxN
        XWX
        Cm % error in coefficients, MxMxD. This is retrieved from the variable cache.
        nonOutlierIndices
        tensionValue
    end

    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Initialization
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function self = MonotonicSmoothingSpline(t,x,distribution,varargin)
            self@SmoothingSpline(t,x,distribution,varargin{:});
        end
    

        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Stuff
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function self = tensionParameterDidChange(self)
            % Tension parameter changed, so we need to recompute the
            % solution, then then compute the PP coefficients for that
            % solution.
            if isa(self.distribution,'NormalDistribution')
                [self.m,~,~,self.isConstrained,self.variableCache] = MonotonicSmoothingSpline.TensionSolution(self.lambda,self.mu,self.t,self.x,self.t_knot,self.K,self.T,self.distribution,1/(self.distribution.sigma^2),self.variableCache);
            elseif ~isempty(self.distribution.w)
                [self.m,~,~,self.isConstrained,self.variableCache] = MonotonicSmoothingSpline.IteratedLeastSquaresTensionSolution(self.lambda,self.mu,self.t,self.x,self.t_knot,self.K,self.T,self.distribution,self.variableCache);
            else
                error('No weight function given! Unable to proceed.');
            end
            
            [self.C,self.t_pp,self.B] = BSpline.PPCoefficientsFromSplineCoefficients( self.m, self.t_knot, self.K, self.B );
            
            self.outlierIndices = find(abs(self.epsilon) > self.outlierThreshold);
        end
        
    end
    methods (Static)
                
        function [m,Cm,CmInv,isConstrained,cachedVars] = TensionSolution(lambda,mu,t,x,t_knot,K,T,distribution,W,cachedVars)
            % N     # of observations
            % M     # of splines
            % Q     # of points in quadrature grid
            %
            % inputs:
            % X         splines on the observation grid, NxM
            % V         spline derivatives on the quadrature grid, QxM
            % sigma     errors of observations, either a scalar, Nx1, OR if 
            %           size(sigma)=[N N], then we assume it's the weight
            %           matrix W (essentially ~diag(1./sigma^2)
            % lambda    tension parameter
            % x         observations (Nx1)
            % mu        mean tension
            % XWX       (optional) precomputed matrix X'*Wx*X
            % XWx       (optional) precomputed matrix X'*Wx*x
            % VV       (optional) precomputed matrix V'*V
            %
            % output:
            % m         coefficients of the splines, Mx1
            % CmInv     Inverse of the covariance of coefficients, MxM
            
            cachedVars = SmoothingSpline.PrecomputeTensionSolutionMatrices(t,x,t_knot,K,T,distribution,W,cachedVars);
            
            N = length(x);
            Q = size(cachedVars.V,1);
            
            XWX = cachedVars.XWX;
            XWx = cachedVars.XWx;
            VV = cachedVars.VV;
 
            E_x = XWX + (lambda*N/Q)*(VV);
            
            % add the mean tension value
            if mu ~= 0.0
                V = cachedVars.V;
                B = XWx + (lambda*N/Q)*mu*transpose(sum( V,1));
            else
                B = XWx;
            end
            
            % Check if the tension is sufficiently high that we just need
            % to use the constrained solution. This scenario only makes
            % sense because we know we've used the interpolation points as
            % knot points... otherwise rcond->0 could mean something else.
            if rcond(E_x) < 5e-14
                if isempty(cachedVars.m_constrained) || isempty(cachedVars.Cm_constrained)
                    t_knot_constrained = cat(1,min(t)*ones(T,1),max(t)*ones(T,1));
                    % T=2 indicates tension on acceleration, which we want
                    % to be zero, so we would want K=2
                    cspline = ConstrainedSpline(t,x,T,t_knot_constrained,distribution,[]);
                    cachedVars.m_constrained = cachedVars.X\cspline(t);
                    
                    % if m_x = inv(X)*Y*m_y then,
                    % Cm_x = G*Cm_y*G^T if G = inv(X)*Y
                    Y = cspline.X;
                    G = cachedVars.X\Y;
                    cachedVars.Cm_constrained = G*(cspline.CmInv\(G.'));
                end
                
                isConstrained = 1;
                m = cachedVars.m_constrained;
                Cm = cachedVars.Cm_constrained;
                CmInv = [];
            else
                isConstrained = 0;
                m = E_x\B;
                Cm = [];
                CmInv = E_x;
            end
            
            cachedVars.Cm = Cm;
            cachedVars.CmInv = CmInv;
        end
        
        % e
        function [m,Cm,CmInv,isConstrained,cachedVars] = IteratedLeastSquaresTensionSolution(lambda,mu,t,x,t_knot,K,T,distribution,cachedVars)
            % Out first call is basically just a normal fit to a tension
            % spline. If W is not set, it will be set to either 1/w0^2 or
            % the correlated version
            cachedVars.W = []; cachedVars.XWX = []; cachedVars.XWx = [];
            [m,Cm,CmInv,isConstrained,cachedVars] = SmoothingSpline.TensionSolution(lambda,mu,t,x,t_knot,K,T,distribution,[],cachedVars);
            
            X = cachedVars.X;
            sigma2_previous = (distribution.sigma0)^2;
            rel_error = 1.0;
            repeats = 1;
            while (rel_error > 0.01)
                sigma_w2 = distribution.w(X*m - x);
                
                if exist('rho_t','var')
                    Sigma2 = (sqrt(sigma_w2) * sqrt(sigma_w2).') .* cachedVars.rho_t;
                    W = inv(Sigma2);
                else
                    W = diag(1./sigma_w2);
                end
                
                % hose any cached variable associated with W...
                cachedVars.W = []; cachedVars.XWX = []; cachedVars.XWx = [];
                % ...and recompute the solution with this new weighting
                [m,Cm,CmInv,isConstrained,cachedVars] = SmoothingSpline.TensionSolution(lambda,mu,t,x,t_knot,K,T,distribution,W,cachedVars);
                
                rel_error = max( abs((sigma_w2-sigma2_previous)./sigma_w2) );
                sigma2_previous=sigma_w2;
                repeats = repeats+1;
                
                if (repeats == 250)
                    disp('Failed to converge after 250 iterations.');
                    break;
                end
            end
            
        end
         
        
        
    end
end
