classdef TensionSpline < BSpline
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        T           % degree at which tension is applied
        lambda      % tension parameter
        mu          % mean value of tension
        w           % weight function
        Cm          % error in coefficients, MxMxD
        
        X
        V
        
        x
        t
        sigma
    end
    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Initialization
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function self = TensionSpline(t,x,lambda,varargin)
            N = length(t);      
            t = reshape(t,[],1); % Nx1
            if size(x,2) == N
                x = x.';
            end
            if size(x,1) ~= N
                error('x and t must have the same length.');
            end
            D = size(x,2);
            
            nargs = length(varargin);
            if mod(nargs,2) ~= 0
                error('Arguments must be given as name/value pairs');
            end
            
            K = 4; % default spline order (cubic spline)
            T = 2; % default tension *degree* (order-1)
            sigma = 1; % default position error
            mu = 0;
            didSetWeightFunction = 0;
            
            for k = 1:2:length(varargin)
                if strcmp(varargin{k}, 'K')
                    K = varargin{k+1};
                elseif strcmp(varargin{k}, 'S')
                    K = varargin{k+1}+1;
                elseif strcmp(varargin{k}, 'T')
                    T = varargin{k+1};
                elseif strcmp(varargin{k}, 'sigma')
                    sigma = varargin{k+1};
                elseif strcmp(varargin{k}, 'mu')
                    mu = varargin{k+1};
                elseif strcmp(varargin{k}, 'weight_function')
                    w = varargin{k+1};
                    didSetWeightFunction = 1;
                end
            end
            
            % Compute the spline values at the observation points
            t_knot = BSpline.KnotPointsForPoints(t,K);
            X = BSpline.Spline( t, t_knot, K, 0 ); % NxM
            
            % Now we need a quadrature (integration) grid that is finer
            Q = 10*N; % number of points on the quadrature grid
            tq = linspace(t(1),t(end),Q)';
            B = BSpline.Spline( tq, t_knot, K, T );
            V = squeeze(B(:,:,T+1)); % QxM
            
            % Now compute the coefficients
            M = size(X,2);
            m = zeros(M,D);
            Cm = zeros(M,M,D);
            for i=1:D
                if didSetWeightFunction == 1
                    [m(:,i),Cm(:,:,i)] = TensionSpline.IteratedLeastSquaresTensionSolution(X,V,sigma,lambda,x(:,i),mu,w);
                else
                    [m(:,i),Cm(:,:,i)] = TensionSpline.TensionSolution(X,V,sigma,lambda,x(:,i),mu);
                end
            end
            
            self@BSpline(t,x,K,t_knot,m);
            self.lambda = lambda;
            self.mu = mu;
            self.T = T;
            self.Cm = Cm;
            self.X = X;
            self.V = V;
            self.t = t;
            self.x = x;
            self.sigma = sigma;
            if didSetWeightFunction == 1
                self.w = w;
            end
            
        end
        
        function B = Splines(self,t,varargin)
            % return the splines being used (evaluated at points t)
            if isempty(varargin) || (length(varargin) == 1 && varargin{1} == 0)
                B = BSpline.Spline( t, self.t_knot, self.K, 0 );
            else
                B = BSpline.Spline( t, self.t_knot, self.K, varargin{1} );
                B = B(:,:,varargin{1});
            end
        end
        
        function self = set.lambda(self,newlambda)
            if self.lambda ~= newlambda
                self.lambda = newlambda;
                self.tensionParameterDidChange();
            end
        end
        
        function self.tensionParameterDidChange(self)
            for i=1:self.D
                if ~isempty(self.w)
                    [self.m(:,i),self.Cm(:,:,i)] = TensionSpline.IteratedLeastSquaresTensionSolution(self.X,self.V,self.sigma,self.lambda,self.x(:,i),self.mu,self.w);
                else
                    [self.m(:,i),self.Cm(:,:,i)] = TensionSpline.TensionSolution(self.X,self.V,self.sigma,self.lambda,self.x(:,i),self.mu);
                end
            end
        end
        
    end
    
    methods (Static)
        function [m, Cm] = TensionSolution(X,V,sigma,lambda,x,mu)
            % N     # of observations
            % M     # of splines
            % Q     # of points in quadrature grid
            %
            % inputs:
            % X         splines on the observation grid, NxM
            % V         spline derivatives on the quadrature grid, QxM
            % sigma     errors of observations, either a scalar, Nx1, OR if size(sigma)=[N N], then we assume it's the weight matrix
            % lambda    tension parameter
            % x         observations (Nx1)
            % mu        mean tension
            %
            % output:
            % m         coefficients of the splines, Mx1
            % Cm        covariance of coefficients, MxM
            N = length(x);
            Q = size(V,1);
            
            if size(sigma,1) == N && size(sigma,2) == N
                Wx = sigma;
                E_x = X'*Wx*X + (lambda*N/Q)*(V'*V); % MxM
                B = X'*Wx*x; % (MxN * NxN * Nx1) = Mx1
            elseif length(sigma) == 1
                E_x = X'*X/(sigma*sigma) + (lambda*N/Q)*(V'*V);  % MxM
                B = X'*x/(sigma*sigma);
            elseif length(sigma) == N
                Wx = diag(1./(sigma.^2));
                E_x = X'*Wx*X + (lambda*N/Q)*(V'*V); % MxM
                B = X'*Wx*x; % (MxN * NxN * Nx1) = Mx1
            else
                error('sigma must have the same length as x and t.');
            end
            
            % add the mean tension value
            if mu ~= 0.0
                B = B + (lambda*N/Q)*mu*transpose(sum( V,1));
            end
            
            % Now s
            m = E_x\B;
            
            if nargout == 2
                Cm = inv(E_x);
            end
        end
        
        function [m,Cm] = IteratedLeastSquaresTensionSolution(X,V,sigma,lambda,x,mu,w)
            if length(sigma) == 1
                sigma = ones(size(x))*sigma;
            end
            Wx = diag(1./(sigma.^2));
            m = TensionSpline.TensionSolution(X,V,Wx,lambda,x,mu);
            
            error_x_previous = sigma*sigma;
            rel_error = 1.0;
            repeats = 1;
            while (rel_error > 0.01)
                dx2 = w(X*m - x);
                
                Wx = diag(1./(dx2));
                
                m = TensionSpline.TensionSolution(X,V,Wx,lambda,x,mu);
                
                rel_error = max( (dx2-error_x_previous)./dx2 );
                error_x_previous=dx2;
                repeats = repeats+1;
                
                if (repeats == 100)
                    disp('Failed to converge after 100 iterations.');
                    break;
                end
            end
            
            
            if nargout == 2
                Cm = inv(X'*Wx*X + (lambda*N/Q)*(V'*V));
            end
        end
    end
    
end

