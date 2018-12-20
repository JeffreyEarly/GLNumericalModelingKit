classdef RobustTensionSpline < TensionSpline
    %RobustTensionSpline Allows for outliers
    
    properties
        indicesOfOutliers = []
        goodIndices = []
        t
        x
        X
        W
    end
    
    properties (Dependent)
       t_all
       x_all
    end
    
    methods
        function self = RobustTensionSpline(t,x,sigma,varargin)
            self@TensionSpline(t,x,sigma,varargin{:});
            self.goodIndices = reshape(1:length(t),[],1);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Smoothing Matrix and Covariance matrix
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function t = get.t(self)
            t = self.t(self.goodIndices);
        end
        
        function x = get.x(self)
            x = self.x(self.goodIndices);
        end
        
        function X = get.X(self)
            X = self.X(self.goodIndices,:);
        end
        
        function W = get.W(self)
            W = self.W(self.goodIndices,self.goodIndices);
        end
        
        function S = SmoothingMatrix(self)
            % The smoothing matrix S takes the observations and maps them
            % onto the estimated true values.
            S = (self.X*self.Cm*self.X.')*self.W;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Measures of error and effective sample size (n_eff)
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function MSE = ExpectedMeanSquareError(self)
            % This is the *expected* mean-square error normalized by the
            % variance. Note that it is a combination of the sample
            % variance and the variance of the mean.
            %
            % From Craven and Wahba, 1979
            
            S = self.SmoothingMatrix;
            SI = (S-eye(size(S)));
            
            MSE = mean((SI*self.x(self.goodIndices)).^2)/(self.sigma*self.sigma) + 2*trace(S)/length(S) - 1;
        end
        
    end
    
    methods (Static)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Methods for solving the least-squares problem
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [m,Cm,W_g,indicesOfOutliers] = IteratedLeastSquaresTensionSolutionWithOutliers(X,V,sigma,lambda,x,mu,w,XWX,XWx,VV)
            % Same calling sequence as the TensionSolution function, but
            % also includes the weight factor, w
            if length(sigma) == 1
                sigma = ones(size(x))*sigma;
            end
            W = diag(1./(sigma.^2));
            m = TensionSpline.TensionSolution(X,V,W,lambda,x,mu,XWX,XWx,VV);
            
            error_x_previous = sigma.*sigma;
            rel_error = 1.0;
            repeats = 1;
            while (rel_error > 0.01)
                dx2 = w(X*m - x);
                
                indicesOfOutliers = find(dx2>1e5*sigma.*sigma);
                goodIndices = setdiff(1:length(x),indicesOfOutliers);
                
                x_g = x(goodIndices);
                dx2_g = dx2(goodIndices);
                X_g = X(goodIndices,:);
                
                W_g = diag(1./(dx2_g));
                
                m = TensionSpline.TensionSolution(X_g,V,W_g,lambda,x_g,mu);
                
                % dropping an outlier counts as 100% error
                if size(dx2_g) == size(error_x_previous)
                    rel_error = max( (dx2_g-error_x_previous)./dx2_g );
                else
                    rel_error = 1;
                end
                error_x_previous=dx2_g;
                repeats = repeats+1;
                
                if (repeats == 100)
                    disp('Failed to converge after 100 iterations.');
                    break;
                end
            end
            
            
            if nargout >= 2
                N = length(x);
                Q = size(V,1);
                Cm = inv(X'*W_g*X + (lambda*N/Q)*(V'*V));
            end
        end
    end
end

