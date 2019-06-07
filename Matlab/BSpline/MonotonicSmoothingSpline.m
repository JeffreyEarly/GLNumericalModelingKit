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

    end
    

    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Initialization
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function self = MonotonicSmoothingSpline(t,x,distribution,varargin)
            self@SmoothingSpline(t,x,distribution,varargin{:});
            self.tensionParameterDidChange();
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
            
            [m0,Cm,CmInv,isConstrained,cachedVars] = SmoothingSpline.TensionSolution(lambda,mu,t,x,t_knot,K,T,distribution,W,cachedVars);
            
            X = cachedVars.X;
            V = cachedVars.V;
            
            N = length(x); % number of points
            M = size(V,2); % number of splines
            Q = size(cachedVars.V,1); % number of quadrature points
             
                        
            % monotonically increasing
%             S = tril(ones(M)); % lower triangle
            
            % monotonically decreasing
            S = -tril(ones(M)); % upper triangle
            S(:,1)=1;
            
            % force existing coefficients to be monotonic
            
            
            xi = optimvar('xi',M,'LowerBound',0,'UpperBound',Inf);
            xi(1).LowerBound = -Inf;
            xi(1).UpperBound = Inf;
            
            XS = X*S;
            VS = V*S;
            
%             f = sum(x.^2) - 2*sum(x.*(XS*xi)) + sum((XS*xi).^2) + (lambda*N/Q/W)*sum((VS*xi).^2);
            f = W*sum( (XS*xi - x).^2 )/N + (lambda/Q)*sum((VS*xi).^2);
            xi0=S\m0;
            
            qprob = optimproblem('Objective',f);
%             opts = optimoptions('quadprog','Algorithm','trust-region-reflective','Display','off');
%             [sol,qfval,qexitflag,qoutput] = solve(qprob,struct('xi',xi0),'options',opts);

opts = optimoptions('quadprog','Algorithm','interior-point-convex','Display','off');
            [sol,qfval,qexitflag,qoutput] = solve(qprob,'options',opts);

            m = S*sol.xi;
            
            XWX = cachedVars.XWX;
            XWx = cachedVars.XWx;
            VV = cachedVars.VV;
 
            E_x = XWX + (lambda*N/Q)*(VV);
            cachedVars.Cm = [];
            cachedVars.CmInv = E_x;
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
