classdef TensionSpline < BSpline
    %TensionSpline Fit noisy data with a tensioned interpolating spline
    %   3 argument initialization
    %       f = TensionSpline(t,x,sigma);
    %   where
    %       t       array of values for the independent axis
    %       x       array of values for the dependent axis 
    %       sigma   std deviation of noise
    %       f       cubic spline interpolant
    %
    %   TensionSpline takes a number of optional input argument pairs.
    %
    %   'weightFunction' the weight function should be specified when the
    %   errors are non-guassian. The weight function should take a single
    %   argument, the error, and map it to the Gaussian equivalent variance
    %   for the error distribution of choice.
    %
    %   'lambda' lambda is the tension parameter, and can be given directly
    %   as a numeric value, or can be a function that takes this
    %   TensionSpline object as an argument, and returns a numeric value.
    %
    % 
    
    properties
        x
        t
        distribution
        
        T           % degree at which tension is applied
        lambda      % tension parameter
        mu          % mean value of tension
        Cm          % error in coefficients, MxMxD
        knot_dof    % knot dofs
        
        X
        V
        W,XWX,XWx,VV

        indicesOfOutliers = []
        goodIndices = []
    end
    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Initialization
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function self = TensionSpline(t,x,distribution,varargin)
            N = length(t);
            t = reshape(t,[],1);
            x = reshape(x,[],1);
            
            if length(x) ~= N
                error('x and t must have the same length.');
            end
            
            if ~isa(distribution,'Distribution')
               error('The distribution must be a Distribution subclass.') 
            end
            
            nargs = length(varargin);
            if mod(nargs,2) ~= 0
                error('Arguments must be given as name/value pairs');
            end
            
            K = 4; % default spline order (cubic spline)
            T = []; % default tension *degree* (order-1)
            mu = 0;
            knot_dof = 1;
            shouldSetKnotDOFAutomatically = 0;
            lambdaArgument = Lambda.optimalIterated;
            indicesOfOutliers = [];
            t_knot = [];
            
            for k = 1:2:length(varargin)
                if strcmp(varargin{k}, 'K')
                    K = varargin{k+1};
                elseif strcmp(varargin{k}, 'S')
                    K = varargin{k+1}+1;
                elseif strcmp(varargin{k}, 'T')
                    T = varargin{k+1};
                elseif strcmp(varargin{k}, 'lambda')
                    lambdaArgument = varargin{k+1};
                elseif strcmp(varargin{k}, 'mu')
                    mu = varargin{k+1};
                elseif strcmp(varargin{k}, 't_knot')
                    t_knot = varargin{k+1};
                elseif strcmp(varargin{k}, 'knot_dof')
                    if ischar(varargin{k+1}) && strcmp(varargin{k+1}, 'auto')
                        shouldSetKnotDOFAutomatically = 1;
                    elseif isnumeric(varargin{k+1}) && varargin{k+1} >= 1
                        knot_dof = varargin{k+1};
                    else
                        error('invalid option for knot_dof. Set to a value >= 1 or auto.');
                    end
                end
            end
            
            if isempty(T)
                T = K-1;
            end
            
            n_eff = [];
            if isenum(lambdaArgument)
                switch lambdaArgument
                    case {Lambda.optimalExpected}
                        u_rms = TensionSpline.EstimateRMSDerivativeFromSpectrum(t,x,sqrt(distribution.variance),1);
                        n_eff = TensionSpline.EffectiveSampleSizeFromUrms(u_rms, t, sqrt(distribution.variance));
                        a_rms = TensionSpline.EstimateRMSDerivativeFromSpectrum(t,x,sqrt(distribution.variance),T);
                        lambda = (n_eff-1)/(n_eff*a_rms.^2);
                    case {Lambda.fullTensionExpected, Lambda.optimalIterated}
                        % if you're going to optimize, it's best to start
                        % near the full tension solution, rather than
                        % (potentially) near zero
                        a_rms = TensionSpline.EstimateRMSDerivativeFromSpectrum(t,x,sqrt(distribution.variance),T);
                        lambda = 1/a_rms.^2;
                end
            elseif isscalar(lambdaArgument)
                lambda = lambdaArgument;
            else
                error('Invalid choice for lambda. Lambda must be either a scalar or the enum Lambda.');
            end
                        
            if shouldSetKnotDOFAutomatically == 1
                if isempty(n_eff)
                    u_rms = TensionSpline.EstimateRMSDerivativeFromSpectrum(t,x,sqrt(distribution.variance),1);
                    n_eff = TensionSpline.EffectiveSampleSizeFromUrms(u_rms, t, sqrt(distribution.variance));
                end
                % conservative estimate
                knot_dof = max(1,floor(0.5*n_eff));
            end
            
            % Compute the spline values at the observation points
            if isempty(t_knot)
                t_knot = InterpolatingSpline.KnotPointsForPoints(t,K,knot_dof);
            end
            X = BSpline.Spline( t, t_knot, K, 0 ); % NxM
            
            % Now we need a quadrature (integration) grid that is finer
            % if S=T we can optimize this much better because it's all
            % piecewise constant.
            Q = 10*N; % number of points on the quadrature grid
            tq = linspace(t(1),t(end),Q)';
            B = BSpline.Spline( tq, t_knot, K, T );
            V = squeeze(B(:,:,T+1)); % QxM
            
            % Precompute some matrices that might be used again later,
            [XWX,XWx,VV] = TensionSpline.PrecomputeTensionSolutionMatrices(X,V,sqrt(distribution.variance),x);
            
            if isa(distribution,'NormalDistribution')
                [m,CmInv] = TensionSpline.TensionSolution(X,V,distribution.sigma,lambda,x,mu,XWX,XWx,VV);
            elseif ~isempty(distribution.w)
                [m,CmInv,W,XWX,XWx,VV] = TensionSpline.IteratedLeastSquaresTensionSolution(X,V,sqrt(distribution.variance),lambda,x,mu,distribution.w,XWX,XWx,VV,t,distribution.rho);
            else
                error('No weight function given! Unable to proceed.');
            end
                        
            self@BSpline(K,t_knot,m);
            self.lambda = lambda;
            self.mu = mu;
            self.T = T;
            self.Cm = inv(CmInv);
            if exist('W','var')
                self.W = W;
            end
            self.X = X;
            self.V = V;
            self.XWX=XWX;
            self.XWx=XWx;
            self.VV=VV;
            self.t = t;
            self.x = x;
            self.knot_dof = knot_dof;
            self.distribution = distribution;
            self.indicesOfOutliers = indicesOfOutliers;
            self.goodIndices = setdiff(1:length(self.x),self.indicesOfOutliers);
            
            if lambdaArgument == Lambda.optimalIterated
                self.minimizeExpectedMeanSquareError();
            end
        end
   
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Stuff
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function B = Splines(self,t,varargin)
            % return the splines being used (evaluated at points t)
            if isempty(varargin) || (length(varargin) == 1 && varargin{1} == 0)
                B = BSpline.Spline( t, self.t_knot, self.K, 0 );
            else
                B = BSpline.Spline( t, self.t_knot, self.K, varargin{1} );
                B = B(:,:,varargin{1}+1);
            end
        end
        
        function [x_T, t_T] = uniqueValuesAtHighestDerivative(self)
            t_T = self.t_knot(self.K:1:(end-self.K+1));
            t_T = t_T(1:end-1) + diff(t_T)/2;
            x_T = self.ValueAtPoints(t_T,self.K-1);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Responding to changes in the tension parameter
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
        function set.lambda(self,newlambda)
            if isempty(self.lambda) 
                self.lambda = newlambda;
            elseif self.lambda ~= newlambda
                self.lambda = newlambda;
                self.tensionParameterDidChange();
            end
        end
        
        function self = tensionParameterDidChange(self)
            % Tension parameter changed, so we need to recompute the
            % solution, then then compute the PP coefficients for that
            % solution.
            if isa(self.distribution,'NormalDistribution')
                [self.m,CmInv] = TensionSpline.TensionSolution(self.X,self.V,self.distribution.sigma,self.lambda,self.x,self.mu,self.XWX,self.XWx,self.VV);
            elseif ~isempty(self.distribution.w)
                % do *not* feed in the previous values of self.XWX,self.XWx
                % when the tension parameter changes for an iterated
                % solution.
                [self.m,CmInv,self.W,self.XWX,self.XWx,self.VV] = TensionSpline.IteratedLeastSquaresTensionSolution(self.X,self.V,sqrt(self.distribution.variance),self.lambda,self.x,self.mu,self.distribution.w,[],[],self.VV,self.t,self.distribution.rho);
            else
                error('No weight function given! Unable to proceed.'); 
            end

            self.Cm = inv(CmInv);
            [self.C,self.t_pp] = BSpline.PPCoefficientsFromSplineCoefficients( self.m, self.t_knot, self.K );
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Smoothing Matrix and Covariance matrix
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function S = smoothingMatrix(self)
            % The smoothing matrix S takes the observations and maps them
            % onto the estimated true values.
            if ~isempty(self.W)
                S = (self.X*self.Cm*self.X.')*self.W;
            else
                S = (self.X*self.Cm*self.X.')/self.distribution.variance;
            end
        end
        
        function S = covarianceMatrixAtPointsForDerivative(self,t,numDerivs)
            % Returns the covariance matrix for a given derivative at the
            % requested points.
            J = self.Splines(t,numDerivs);            
            if ~isempty(self.distribution.w)
                S = (J*self.Cm*J.')*self.W;
            else
                S = (J*self.Cm*J.')/self.distribution.variance;
            end
        end    
        
        function S = covarianceMatrixForDerivative(self,numDerivs)
            % Returns the covariance matrix for a given derivative at the
            % points of observation
            S = self.covarianceMatrixAtPointsForDerivative(self.t,numDerivs);
        end
        
        function SE = standardErrorAtPointsForDerivative(self,t,numDerivs)
            % Returns the standard error for a given derivative at the
            % points requested.
            SE = sqrt(diag(self.covarianceMatrixAtPointsForDerivative(t,numDerivs)));
        end
        
        function SE = standardErrorForDerivative(self,numDerivs)
            % Returns the standard error for a given derivative at the
            % observation points.
           SE = sqrt(diag(self.covarianceMatrixForDerivative(numDerivs)));
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Measures of error and effective sample size (n_eff)
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function epsilon = epsilon(self)
            epsilon = self.x - self.ValueAtPoints(self.t);
        end
        
        function X2 = sampleVariance(self)
            % The sample variance. Same as ||(S-I)*x||^2/N
            % Not normalized by the variance.
            X2 = mean( self.epsilon.^2,1);
        end
        
        function X2 = sampleVarianceInRange(self,zmin,zmax)
            epsilon = self.epsilon;
            X2 = mean( epsilon( epsilon >= zmin & epsilon <= zmax ).^2,1);
        end
        
        function X2 = sampleVarianceInPercentileRange(self,pctmin,pctmax)
            zmin = self.distribution.locationOfCDFPercentile(pctmin);
            zmax = self.distribution.locationOfCDFPercentile(pctmax);
            X2 = self.sampleVarianceInRange(zmin,zmax);
        end
        
        function SE2 = varianceOfTheMean(self)
            % The variance of the mean is the square of the standard error
            % Normalized by the variance.
            S = self.smoothingMatrix;
            SE2 = self.distribution.variance*trace(S)/length(S);
        end
        
        function MSE = expectedMeanSquareError(self)
            % This is the *expected* mean-square error normalized by the
            % variance. Note that it is a combination of the sample
            % variance and the variance of the mean.
            %
            % From Craven and Wahba, 1979
            
            if ~isempty(self.XWX)
                % S = X*C*X'*W, so trace(S)=trace(C*X'*W*X) under the
                % cyclic property
                MSE = self.sampleVariance/self.distribution.variance + 2*trace(squeeze(self.Cm(:,:,1))*self.XWX)/length(self.x) - 1;
            else
                S = self.smoothingMatrix;
                SI = (S-eye(size(S)));
                
                MSE = mean((SI*self.x).^2)/self.distribution.variance + 2*trace(S)/length(S) - 1;
            end
        end
        
        function [MSE, n] = expectedMeanSquareErrorInRange(self,zmin,zmax)
            epsilon = self.epsilon;
            indices = find(epsilon >= zmin & epsilon <= zmax);
            X2 = mean(epsilon(indices).^2,1);
            expectedVariance = self.distribution.varianceInRange(zmin,zmax);
            
            S = self.smoothingMatrix;
            S = S(indices,indices); 
            
            n = length(S);
            MSE = X2/expectedVariance + 2*trace(S)/n - 1;
        end
        
        function MSE = expectedMeanSquareErrorFromGCV(self)
            % generalized cross-validation estimate for the mean square
            % error from Craven and Wahba, 1979 equation 1.9.
            S = self.smoothingMatrix;
            SI = (S-eye(size(S)));
            a = mean((SI*self.x).^2);
            b = trace(S)/length(S);
            
            MSE = a/(1-b)^2;
            
%             MSE = a*(b/(1-b));
            %             
            %             fprintf(' sigma^2=%.2f ',a/(1-b));
        end
        
        % This MSE is slightly higher than what we actually get, increase
        % as a function of derivative.
        function [MSE,noise] = ExpectedMeanSquareErrorAtAllOrders(self)
           MSE = zeros(self.K,1);
           noise = zeros(self.K,1);
           
           S = self.smoothingMatrix;
           SI = (S-eye(size(S)));
           sigma2 = self.distribution.variance;
           MSE(1) = mean((SI*self.x).^2) + sigma2*(2*trace(S)/length(S) - 1);
           noise(1) = self.distribution.variance;
           
           for iDiff=1:(self.K-1)
               Diff = TensionSpline.FiniteDifferenceMatrixNoBoundary(iDiff,self.t,1);
               DS = Diff*S;
               DminusDS = Diff - DS;
               MSE(iDiff+1) = (sum((DminusDS*self.x).^2) - sigma2*sum(sum(DminusDS.^2)) + sigma2*sum(sum(DS.^2)))/length(Diff);
%                MSE(iDiff+1) = ( sum((DminusDS*self.x).^2) + 2*sigma2*sum(sum(Diff.*DS)) - sigma2*sum(sum(Diff.^2)) )/length(Diff);
               noise(iDiff+1) = sigma2*sum(sum(Diff.^2))/length(Diff);
               % This is the equivalent of sigma2 for the derivative of the noise
           end
        end
        
        function SNR = SignalToNoiseRatio(self)
            sigma2 = self.distribution.variance;
            x_smoothed = self.ValueAtPoints(self.t);
            
            SNR = zeros(self.K,1);
            SNR(1) = mean(x_smoothed.^2)/sigma2;
            for iDiff=1:(self.K-1)
               Diff = TensionSpline.FiniteDifferenceMatrixNoBoundary(iDiff,self.t,1);
               SNR(1+iDiff) = sum((Diff*x_smoothed).^2)/(sigma2*sum(sum(Diff.^2)) );
            end
        end
        
        function SSER = SignalToStandardErrorRatio(self)
            sigma2 = self.distribution.variance;
            x_smoothed = self.ValueAtPoints(self.t);
            S = self.smoothingMatrix;
            
            SSER = zeros(self.K,1);
            SSER(1) = mean(x_smoothed.^2)/sigma2;
            for iDiff=1:(self.K-1)
               Diff = TensionSpline.FiniteDifferenceMatrixNoBoundary(iDiff,self.t,1);
               SE = sigma2*sum(sum( Diff.*(Diff*S) ));
               SSER(1+iDiff) = sum((Diff*x_smoothed).^2)/SE;
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Measures of effective sample size ESS
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function n_eff = effectiveSampleSizeFromVarianceOfTheMean(self)
            n_eff = self.distribution.variance ./ self.varianceOfTheMean;
        end
        
        function n_eff = effectiveSampleSizeFromSampleVariance(self)
            n_eff = 1./(1-self.sampleVariance/self.distribution.variance);
        end

        function n_eff = effectiveSampleSizeFromExpectedMeanSquareError(self)
            n_eff = 1./self.expectedMeanSquareError;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Likelihood
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function phi = LogLikelihood(self)
            Q = 10*length(self.t); % number of points on the quadrature grid
            tq = linspace(self.t(1),self.t(end),Q)';
            phi = self.sampleVariance/(self.sigma^2) + self.lambda *  mean( self.ValueAtPoints(tq,self.T).^2 );
        end
        
        function AIC = AIC(self)
            m = size(self.Cm,1);
            k = self.K * m;
            AIC = 2 * k - 2*self.LogLikelihood();
        end
        
        function AICc = AICc(self)
            m = size(self.Cm,1);
            k = self.K * m;
            n = length(self.x);
            AICc = self.AIC + (2*k*k + 2*k)/(n-k-1);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Minimization
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function lambda = minimize(self,penaltyFunction)
            % the penalty function *must* take a tension spline object and
            % return a scalar. The function will be minimized by varying
            % lambda.
            lambda = TensionSpline.minimizeFunctionOfSpline(self,penaltyFunction);
        end
        
        function lambda = minimizeExpectedMeanSquareError(self)
            lambda = self.minimize(@(aTensionSpline) aTensionSpline.expectedMeanSquareError);
        end
        
        function lambda = minimizeExpectedMeanSquareErrorInPercentileRange(self,pctmin,pctmax)
            zmin = self.distribution.locationOfCDFPercentile(pctmin);
            zmax = self.distribution.locationOfCDFPercentile(pctmax); 
            lambda = self.minimize( @(aTensionSpline) aTensionSpline.expectedMeanSquareErrorInRange(zmin,zmax) );
        end
        
        function lambda = minimizeMeanSquareError(self,t_true,x_true)
           mse = @(aTensionSpline) mean((aTensionSpline.ValueAtPoints(t_true)-x_true).^2)/(aTensionSpline.distribution.variance);
           lambda = self.minimize(mse);
        end
    end
    
    methods (Static)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Minimization
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % We wrap and then unwrap the user-provided penalty function so
        % that it take log10(lambda) rather than lambda. This is much
        % better for the fminsearch algorithm.

        function lambda = minimizeFunctionOfSpline(aTensionSpline,functionOfSpline)
            epsilon = 1e-15;
            errorFunction = @(log10lambdaPlusEpsilon) TensionSpline.FunctionOfSplineWrapper(aTensionSpline,log10lambdaPlusEpsilon,functionOfSpline);
            optimalLog10lambdaPlusEpsilon = fminsearch( errorFunction, log10(aTensionSpline.lambda+epsilon), optimset('TolX', 0.01, 'TolFun', 0.01) );
            lambda = 10^optimalLog10lambdaPlusEpsilon - epsilon;
        end
        
        function error = FunctionOfSplineWrapper(aTensionSpline, log10lambdaPlusEpsilon, functionOfSpline)
            epsilon = 1e-15;
            aTensionSpline.lambda = 10^log10lambdaPlusEpsilon-epsilon;              
            error = functionOfSpline(aTensionSpline);
        end
        
        % Same thing as above, but for a bunch of splines given as a cell
        % array
        function lambda = minimizeFunctionOfSplines(tensionSplines,functionOfSplines)
            epsilon = 1e-15;
            errorFunction = @(log10lambdaPlusEpsilon) TensionSpline.functionOfSplinesWrapper(tensionSplines,log10lambdaPlusEpsilon,functionOfSplines);
            optimalLog10lambdaPlusEpsilon = fminsearch( errorFunction, log10(tensionSplines{1}.lambda+epsilon), optimset('TolX', 0.01, 'TolFun', 0.01) );
            lambda = 10^optimalLog10lambdaPlusEpsilon - epsilon;
        end
        
        function error = functionOfSplinesWrapper(tensionSplines, log10lambdaPlusEpsilon, functionOfSplines)
            epsilon = 1e-15;
            for iSpline=1:length(tensionSplines)
                tensionSplines{iSpline}.lambda = 10^log10lambdaPlusEpsilon-epsilon;
            end
            error = functionOfSplines(tensionSplines);
        end
        
        function error = expectedMeanSquareErrorOfSplines(tensionSplines,zmin,zmax)
            MSE = 0;
            N = 0;
            for iSpline=1:length(tensionSplines)
                [iMSE, n] = tensionSplines{iSpline}.expectedMeanSquareErrorInRange(zmin,zmax);
                MSE = MSE + n*iMSE;
                N = N + n;
            end
            error = MSE/N;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Methods for solving the least-squares problem
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [XWX,XWx,VV] = PrecomputeTensionSolutionMatrices(X,V,sigma,x)
            N = length(x);
            if size(sigma,1) == N && size(sigma,2) == N
                XWX = X'*sigma*X;
                XWx = X'*sigma*x;
            elseif length(sigma) == 1
                XWX = X'*X/(sigma*sigma);
                XWx = X'*x/(sigma*sigma);
            elseif length(sigma) == N
                XWX = X'*diag(1./(sigma.^2))*X; % MxM
                XWx = X'*diag(1./(sigma.^2))*x; % (MxN * NxN * Nx1) = Mx1
            else
                error('sigma must have the same length as x and t.');
            end
            VV = V'*V;
        end
        
        function [m,CmInv,XWX,XWx,VV] = TensionSolution(X,V,sigma,lambda,x,mu,XWX,XWx,VV)
            % N     # of observations
            % M     # of splines
            % Q     # of points in quadrature grid
            %
            % inputs:
            % X         splines on the observation grid, NxM
            % V         spline derivatives on the quadrature grid, QxM
            % sigma     errors of observations, either a scalar, Nx1, OR if 
            %           size(sigma)=[N N], then we assume it's the weight
            %           matrix W
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
            N = length(x);
            Q = size(V,1);
            
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
            
            if ~exist('VV','var') || isempty(VV)
                VV = V'*V;
            end
            
            E_x = XWX + (lambda*N/Q)*(VV);
                        
            % add the mean tension value
            if mu ~= 0.0
                B = XWx + (lambda*N/Q)*mu*transpose(sum( V,1));
            else
                B = XWx;
            end
            
            % Now solve
            m = E_x\B;
            
            CmInv = E_x;
        end
        
        function [m,CmInv,W,XWX,XWx,VV] = IteratedLeastSquaresTensionSolution(X,V,sigma,lambda,x,mu,w,XWX,XWx,VV,t,rho)
            % Same calling sequence as the TensionSolution function, but
            % also includes the weight factor, w
            if exist('t','var') && exist('rho','var') && ~isempty(rho)
                rho_t = rho(t - t.');   
                if length(sigma) == 1
                    sigma = ones(size(x))*sigma;
                end      
                Sigma2 = (sigma * sigma.') .* rho_t;
                W = inv(Sigma2);
                m = TensionSpline.TensionSolution(X,V,W,lambda,x,mu,XWX,XWx,VV);
            else
                m = TensionSpline.TensionSolution(X,V,sigma,lambda,x,mu,XWX,XWx,VV);
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
                
                [m,CmInv,XWX,XWx,VV] = TensionSpline.TensionSolution(X,V,W,lambda,x,mu);
                
                rel_error = max( (sigma_w2-sigma2_previous)./sigma_w2 );
                sigma2_previous=sigma_w2;
                repeats = repeats+1;
                
                if (repeats == 100)
                    disp('Failed to converge after 100 iterations.');
                    break;
                end
            end
            
        end
        
        
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
                N = length(x_g);
                Q = size(V,1);
                Cm = inv(X_g'*W_g*X_g + (lambda*N/Q)*(V'*V));
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Supporting function
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [sigma, mu] = StandardDeviationAndMeanOfInterquartileRange(x)
            mu = median(x);
            x = sort(x-mu);
            x_Q1 = median(x(1:floor(length(x)/2)));
            x_Q3 = median(x(ceil(length(x)/2)+1:end));
            sigma = (x_Q3-x_Q1)/1.349;
        end
        
        function sigma = StandardDeviationOfInterquartileRange(x)
            x = sort(x);
            x_Q1 = median(x(1:floor(length(x)/2)));
            x_Q3 = median(x(ceil(length(x)/2)+1:end));
            sigma = (x_Q3-x_Q1)/1.349;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Methods for estimating lambda
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
        function n_eff = EffectiveSampleSizeFromUrms(u_rms, t, sigma)
            % These are the coefficients of the empirical best fits for
            % slopes [-2,-3,-4] to the model n_eff = exp(b)*gamma^m
            m = [0.6652; 0.7904; 0.8339];
            b = [1.6903; 2.1786; 2.3288];
            
            dt = median(abs(diff(t)));
            
            gamma = sigma/( u_rms*dt );
            n_eff = max(1,exp(b(1))*gamma.^m(1));
        end
                        
        function [a_rms, a_std, a_mean] = EstimateRMSDerivativeFromSpectrum( t, x, sigma, D, shouldPlotSpectra)
            % EstimateRMSDerivativeFromSpectrum Given some signal (t,x)
            % contaminated by noise sigma, this uses the spectrum to
            % estimate u_rms.
            %
            % D = 1 is velocity, D=2 is acceleration
            
            xin = x;
            tin = t;
            
            if TensionSpline.IsEvenlySampled(t) ~= 1
                %    fprintf('interpolating...\n');
                dt = round(median(diff(t)));
                N = ceil((t(end)-t(1))/dt);
                t2 = dt*((0:(N-1))') + t(1);
                x = interp1(t,x,t2);
                t = t2;
            end
            
            [p,~,mu]=polyfit(t,x,D);
            a_mean = factorial(D)*p(1)/mu(2)^D;
            
            % now remove the trend
            x = x-polyval(p,t,[],mu);
            
            if 1 == 0
                dt = t(2) - t(1);
                T = t(end)-t(1);
                N = length(t);
                
                df = 1/T;
                f = ([0:ceil(N/2)-1 -floor(N/2):-1]*df)';
                
                ubar = fft(x);
                s_signal = (ubar.*conj(ubar)) .* (2*pi*f).^(2*D) * (dt/N);
            else
                [DiffMatrix,t_u] = TensionSpline.FiniteDifferenceMatrixNoBoundary(D,t,1);
                
                dt = t_u(2)-t_u(1);
                T = t_u(end)-t_u(1);
                N = length(t_u);
                
                df = 1/T;
                f = ([0:ceil(N/2)-1 -floor(N/2):-1]*df)';
                
                ubar = fft(DiffMatrix*x);
                s_signal = (ubar .* conj(ubar)) * (dt/N);

    % This is a one-sided spectrum, need to double the cutoff below.            
%                 psi = sleptap(size(t_u,1),1);
%                 [f,s_signal] = mspec(t_u(2)-t_u(1),DiffMatrix*xin,psi,'cyclic');
%                 df = f(2)-f(1);
            end
            s_noise = sigma*sigma*dt*(2*pi*f).^(2*D);
            
            % The factor of 10 is consistent with 80% confidence.
            % 95% confidence (actually, 97.5% ?) is 39.5
            alpha = 0.10;
            cutoff = 2/TensionSpline.chi2inv(alpha/2,2);
            
            % There are two ways to think of this. Either you look for the
            % 95% range of the signal, or the 95% range of the expected
            % noise.
            alpha = 0.99999;
            dof = 2;
            cutoff = 1*TensionSpline.chi2inv(alpha,dof)/dof;
            
            u2 = sum((s_signal > cutoff*s_noise) .* s_signal)*df;
            a_std = sqrt(u2);
            a_rms = sqrt( u2 + a_mean^2 );
            
            if nargin > 4 && shouldPlotSpectra == 1
                f = fftshift(f);
                s_signal = fftshift(s_signal);
                s_noise = fftshift(s_noise);
                
                figure
                plot(f,s_signal)
                hold on
                plot(f,cutoff*s_noise), ylog
                
%                 figure
%                 plot(tin,xin), hold on, plot(tin,polyval(p,tin,[],mu))
            end
        end
        
        function flag = IsEvenlySampled(t)
            % Checks the sampling rate of t.
            % Returns 1 if the data is evenly sampled (a single unique dt)
            % Returns 2 if the data is sampled with multiples of a unique dt
            % Return 0 otherwise
            unique_dt = unique(diff(t));
            if length(unique_dt) == 1
                flag = 1;
            else
                dt_multiples = unique_dt/min(unique_dt);
                if all( dt_multiples-1.0 < 1e-7 )
                    flag = 1;
                elseif all(mod(dt_multiples,1.0) < 0.01)
                    flag = 2;
                else
                    flag = 0;
                end
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Optimizing lambda when the true values are known
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        function [rmse,norm] = MeanSquareErrorAtAllOrders(aTensionSpline, t_true, x_true)
            if length(unique(diff(t_true))) > 1
                error('This only works for evenly spaced data at the moment.');
            end
            
            rmse = zeros(aTensionSpline.K,1);
            norm = zeros(aTensionSpline.K,1);
            dt = t_true(2)-t_true(1);
            for D = 0:(aTensionSpline.K-1)
                
                % differentiate D times
                % We remove the mean from *position*
                signal_true = x_true;
                for i=1:D
                    signal_true = diff(signal_true)/dt;
                end
                
                if 1 == 0
                     % points of accuracy move dt/2 further inside each D
                    t_signal = t_true(1:(length(t_true)-D)) + D*dt/2;
                    
                    signal_obs = aTensionSpline.ValueAtPoints(t_signal,D);
                else
                    signal_obs = aTensionSpline.ValueAtPoints(t_true);
                    for i=1:D
                        signal_obs = diff(signal_obs)/dt;
                    end
                end
                
                rmse(D+1) = sqrt( mean( (signal_true - signal_obs).^2 ) );
                norm(D+1) = sqrt(mean(signal_true.^2));
                if D == 0
                    norm(D+1) = sqrt(mean((signal_true-mean(signal_true)).^2));
                end
            end
        end
        
        function Q = QScore(aTensionSpline, t_true, x_true)
            [rmse,norm] = TensionSpline.MeanSquareErrorAtAllOrders(aTensionSpline, t_true, x_true);
            Q = mean(rmse./norm);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Supporting finite difference routines
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        function [D,z,width] = FiniteDifferenceMatrixNoBoundary(numDerivs, x, order)
            % Creates a finite difference matrix of aribtrary accuracy, on an arbitrary
            % grid. It does not implement boundary conditions (check my other routine
            % for that), because it seeks to make all rows linearly independent.
            %
            %
            % numDerivs ? the number of derivatives
            % x ? the grid
            % z location where approximations are to be accurate,
            % orderOfAccuracy ? minimum order of accuracy required
            % width ? the distance between the first and last point used in the
            % approximation.
            %
            % Jeffrey J. Early, 2015
            
            n = length(x);
            m = n - numDerivs;
            D = zeros(m,n);
            width = zeros(m,1);
            
            % order != accurracy.
            nPoints = (numDerivs+1) + 2*(order-1);
            
            if mod(numDerivs,2) == 0
                half = numDerivs/2;
                z = x((1+half):(n-half));
            else
                mids = x(1:(n-1)) + diff(x)/2;
                half = floor(numDerivs/2);
                z = mids((1+half):(end-half));
            end
            
            % do we want to find the n closest points?
            for i=1:m
                
                range_left = find( x <= z(i), ceil(nPoints/2), 'last');
                range_right = find( x > z(i), nPoints - length(range_left), 'first');
                range = union(range_left,range_right);
                
                if length(range)<nPoints
                    range_right = find( x >= z(i), ceil(nPoints/2), 'first');
                    range_left = find( x < z(i), nPoints - length(range_right), 'last');
                    range = union(range_left,range_right);
                end
                
                c = TensionSpline.weights( z(i), x(range), numDerivs );
                D(i,range) = c(numDerivs+1,:);
                width(i) = max(x(range))-min(x(range));
            end
            
        end
        
        function c = weights(z,x,m)
            % Calculates FD weights. The parameters are:
            %  z   location where approximations are to be accurate,
            %  x   vector with x-coordinates for grid points,
            %  m   highest derivative that we want to find weights for
            %  c   array size m+1,lentgh(x) containing (as output) in
            %      successive rows the weights for derivatives 0,1,...,m.
            %
            % Taken from Bengt Fornberg
            %
            n=length(x); c=zeros(m+1,n); c1=1; c4=x(1)-z; c(1,1)=1;
            for i=2:n
                mn=min(i,m+1); c2=1; c5=c4; c4=x(i)-z;
                for j=1:i-1
                    c3=x(i)-x(j);  c2=c2*c3;
                    if j==i-1
                        c(2:mn,i)=c1*((1:mn-1)'.*c(1:mn-1,i-1)-c5*c(2:mn,i-1))/c2;
                        c(1,i)=-c1*c5*c(1,i-1)/c2;
                    end
                    c(2:mn,j)=(c4*c(2:mn,j)-(1:mn-1)'.*c(1:mn-1,j))/c3;
                    c(1,j)=c4*c(1,j)/c3;
                end
                c1=c2;
            end
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Supporting functions for chi2inv
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        function x = chi2inv(p,v)
            %CHI2INV Inverse of the chi-square cumulative distribution function (cdf).
            %   X = CHI2INV(P,V)  returns the inverse of the chi-square cdf with V
            %   degrees of freedom at the values in P. The chi-square cdf with V
            %   degrees of freedom, is the gamma cdf with parameters V/2 and 2.
            %
            %   The size of X is the common size of P and V. A scalar input
            %   functions as a constant matrix of the same size as the other input.
            %
            %   See also CHI2CDF, CHI2PDF, CHI2RND, CHI2STAT, ICDF.
            
            %   References:
            %      [1]  M. Abramowitz and I. A. Stegun, "Handbook of Mathematical
            %      Functions", Government Printing Office, 1964, 26.4.
            %      [2] E. Kreyszig, "Introductory Mathematical Statistics",
            %      John Wiley, 1970, section 10.2 (page 144)
            
            %   Copyright 1993-2004 The MathWorks, Inc.
            %   $Revision: 2.10.2.2 $  $Date: 2004/07/05 17:02:22 $
            
            % Call the gamma inverse function.
            x = TensionSpline.gaminv(p,v/2,2);
            
            % Return NaN if the degrees of freedom is not positive.
            k = (v <= 0);
            if any(k(:))
                x(k) = NaN;
            end
        end
        
        function [x,xlo,xup] = gaminv(p,a,b,pcov,alpha)
            %GAMINV Inverse of the gamma cumulative distribution function (cdf).
            %   X = GAMINV(P,A,B) returns the inverse cdf for the gamma distribution
            %   with shape A and scale B, evaluated at the values in P.  The size of X
            %   is the common size of the input arguments.  A scalar input functions as
            %   a constant matrix of the same size as the other inputs.
            %
            %   [X,XLO,XUP] = GAMINV(P,A,B,PCOV,ALPHA) produces confidence bounds for
            %   X when the input parameters A and B are estimates. PCOV is a 2-by-2
            %   matrix containing the covariance matrix of the estimated parameters.
            %   ALPHA has a default value of 0.05, and specifies 100*(1-ALPHA)%
            %   confidence bounds.  XLO and XUP are arrays of the same size as X
            %   containing the lower and upper confidence bounds.
            %
            %   See also GAMCDF, GAMFIT, GAMLIKE, GAMPDF, GAMRND, GAMSTAT.
            
            %   GAMINV uses Newton's method to find roots of GAMCDF(X,A,B) = P.
            
            %   References:
            %      [1] Abramowitz, M. and Stegun, I.A. (1964) Handbook of Mathematical
            %          Functions, Dover, New York, section 26.1.
            %      [2] Evans, M., Hastings, N., and Peacock, B. (1993) Statistical
            %          Distributions, 2nd ed., Wiley.
            
            %   Copyright 1993-2004 The MathWorks, Inc.
            %   $Revision: 2.10.2.8 $  $Date: 2004/07/28 04:38:22 $
            
            if nargin < 2
                error('stats:gaminv:TooFewInputs',...
                    'Requires at least two input arguments.');
            elseif nargin < 3
                b = 1;
            end
            
            % More checking if we need to compute confidence bounds.
            if nargout > 2
                if nargin < 4
                    error('stats:gaminv:TooFewInputs',...
                        'Must provide covariance matrix to compute confidence bounds.');
                end
                if ~isequal(size(pcov),[2 2])
                    error('stats:gaminv:BadCovariance',...
                        'Covariance matrix must have 2 rows and columns.');
                end
                if nargin < 5
                    alpha = 0.05;
                elseif ~isnumeric(alpha) || numel(alpha) ~= 1 || alpha <= 0 || alpha >= 1
                    error('stats:gaminv:BadAlpha',...
                        'ALPHA must be a scalar between 0 and 1.');
                end
            end
            
            % Weed out any out of range parameters or edge/bad probabilities.
            try
                okAB = (0 < a) & (0 < b);
                k = (okAB & (0 < p & p < 1));
            catch
                error('stats:gaminv:InputSizeMismatch',...
                    'Non-scalar arguments must match in size.');
            end
            allOK = all(k(:));
            
            % Fill in NaNs for out of range cases, fill in edges cases when P is 0 or 1.
            if ~allOK
                x = repmat(NaN, size(k));
                x(okAB & p == 0) = 0;
                x(okAB & p == 1) = Inf;
                
                
                if nargout > 1
                    xlo = x; % NaNs or zeros or Infs
                    xup = x; % NaNs or zeros or Infs
                end
                
                % Remove the bad/edge cases, leaving the easy cases.  If there's
                % nothing remaining, return.
                if any(k(:))
                    if numel(p) > 1, p = p(k); end
                    if numel(a) > 1, a = a(k); end
                    if numel(b) > 1, b = b(k); end
                else
                    return;
                end
            end
            
            % ==== Newton's Method to find a root of GAMCDF(X,A,B) = P ====
            
            % Limit this to maxiter iterations.
            maxiter = 500;
            iter = 0;
            
            % Choose a starting guess for q.  Use quantiles from a lognormal
            % distribution with the same mean (==a) and variance (==a) as G(a,1).
            loga = log(a);
            sigsq = log(1+a) - loga;
            mu = loga - 0.5 .* sigsq;
            q = exp(mu - sqrt(2.*sigsq).*erfcinv(2*p));
            
            h = ones(size(p));
            
            % Break out of the iteration loop when the relative size of the last step
            % is small for all elements of q.
            reltol = eps(class(p)).^(3/4);
            dF = zeros(size(p));
            while any(abs(h(:)) > reltol*q(:))
                iter = iter + 1;
                if iter > maxiter
                    % Too many iterations.  This should not happen.
                    break
                end
                
                F = TensionSpline.gamcdf(q,a,1);
                f = max(TensionSpline.gampdf(q,a,1), realmin(class(p)));
                dF = F-p;
                h = dF ./ f;
                qnew = q - h;
                % Make sure that the current iterates stay positive.  When Newton's
                % Method suggests steps that lead to negative values, take a step
                % 9/10ths of the way to zero instead.
                ksmall = find(qnew <= 0);
                if ~isempty(ksmall)
                    qnew(ksmall) = q(ksmall) / 10;
                    h = q - qnew;
                end
                q = qnew;
            end
            
            badcdf = (isfinite(a(:)) & abs(dF(:))>sqrt(eps));
            if iter>maxiter || any(badcdf)   % too many iterations or cdf is too far off
                didnt = find(abs(h)>reltol*q | badcdf);
                didnt = didnt(1);
                if numel(a) == 1, abad = a; else abad = a(didnt); end
                if numel(b) == 1, bbad = b; else bbad = b(didnt); end
                if numel(p) == 1, pbad = p; else pbad = p(didnt); end
                warning('stats:gaminv:NoConvergence',...
                    'GAMINV did not converge for a = %g, b = %g, p = %g.',...
                    abad,bbad,pbad);
            end
            
            % Add in the scale factor, and broadcast the values to the correct place if
            % need be.
            if allOK
                x = q .* b;
            else
                x(k) = q .* b;
            end
            
            % Compute confidence bounds if requested.
            if nargout >= 2
                logq = log(q);
                dqda = -dgammainc(q,a) ./ exp((a-1).*logq - q - gammaln(a));
                
                % Approximate the variance of x=q*b on the log scale.
                %    dlogx/da = dlogx/dq * dqda = dqda/q
                %    dlogx/db = 1/b
                logx = logq + log(b);
                varlogx = pcov(1,1).*(dqda./q).^2 + 2.*pcov(1,2).*dqda./(b.*q) + pcov(2,2)./(b.^2);
                if any(varlogx(:) < 0)
                    error('stats:gaminv:BadCovariance',...
                        'PCOV must be a positive semi-definite matrix.');
                end
                z = -norminv(alpha/2);
                halfwidth = z * sqrt(varlogx);
                
                % Convert back to original scale
                if allOK
                    xlo = exp(logx - halfwidth);
                    xup = exp(logx + halfwidth);
                else
                    xlo(k) = exp(logx - halfwidth);
                    xup(k) = exp(logx + halfwidth);
                end
            end
        end
        
        function [p,plo,pup] = gamcdf(x,a,b,pcov,alpha)
            %GAMCDF Gamma cumulative distribution function.
            %   P = GAMCDF(X,A,B) returns the gamma cumulative distribution function
            %   with shape and scale parameters A and B, respectively, at the values in
            %   X.  The size of P is the common size of the input arguments.  A scalar
            %   input functions as a constant matrix of the same size as the other
            %   inputs.
            %
            %   Some references refer to the gamma distribution with a single
            %   parameter.  This corresponds to the default of B = 1.
            %
            %   [P,PLO,PUP] = GAMCDF(X,A,B,PCOV,ALPHA) produces confidence bounds for
            %   P when the input parameters A and B are estimates.  PCOV is a 2-by-2
            %   matrix containing the covariance matrix of the estimated parameters.
            %   ALPHA has a default value of 0.05, and specifies 100*(1-ALPHA)%
            %   confidence bounds.  PLO and PUP are arrays of the same size as P
            %   containing the lower and upper confidence bounds.
            %
            %   See also GAMFIT, GAMINV, GAMLIKE, GAMPDF, GAMRND, GAMSTAT.
            
            %   GAMMAINC does computational work.
            
            %   References:
            %      [1] Abramowitz, M. and Stegun, I.A. (1964) Handbook of Mathematical
            %          Functions, Dover, New York, section 26.1.
            %      [2] Evans, M., Hastings, N., and Peacock, B. (1993) Statistical
            %          Distributions, 2nd ed., Wiley.
            
            %   Copyright 1993-2004 The MathWorks, Inc.
            %   $Revision: 2.12.2.3 $  $Date: 2004/01/24 09:33:52 $
            
            if nargin < 2
                error('stats:gamcdf:TooFewInputs',...
                    'Requires at least two input arguments.');
            elseif nargin < 3
                b = 1;
            end
            
            % More checking if we need to compute confidence bounds.
            if nargout > 1
                if nargin < 4
                    error('stats:gamcdf:TooFewInputs',...
                        'Must provide covariance matrix to compute confidence bounds.');
                end
                if ~isequal(size(pcov),[2 2])
                    error('stats:gamcdf:BadCovariance',...
                        'Covariance matrix must have 2 rows and columns.');
                end
                if nargin < 5
                    alpha = 0.05;
                elseif ~isnumeric(alpha) || numel(alpha) ~= 1 || alpha <= 0 || alpha >= 1
                    error('stats:gamcdf:BadAlpha',...
                        'ALPHA must be a scalar between 0 and 1.');
                end
            end
            
            % Return NaN for out of range parameters.
            a(a <= 0) = NaN;
            b(b <= 0) = NaN;
            x(x < 0) = 0;
            
            try
                z = x ./ b;
                p = gammainc(z, a);
            catch
                error('stats:gamcdf:InputSizeMismatch',...
                    'Non-scalar arguments must match in size.');
            end
            p(z == Inf) = 1;
            
            % Compute confidence bounds if requested.
            if nargout >= 2
                % Approximate the variance of p on the logit scale
                logitp = log(p./(1-p));
                dp = 1 ./ (p.*(1-p)); % derivative of logit(p) w.r.t. p
                da = dgammainc(z,a) .* dp; % dlogitp/da = dp/da * dlogitp/dp
                db = -exp(a.*log(z)-z-gammaln(a)-log(b)) .* dp; % dlogitp/db = dp/db * dlogitp/dp
                varLogitp = pcov(1,1).*da.^2 + 2.*pcov(1,2).*da.*db + pcov(2,2).*db.^2;
                if any(varLogitp(:) < 0)
                    error('stats:gamcdf:BadCovariance',...
                        'PCOV must be a positive semi-definite matrix.');
                end
                
                % Use a normal approximation on the logit scale, then transform back to
                % the original CDF scale
                halfwidth = -norminv(alpha/2) * sqrt(varLogitp);
                explogitplo = exp(logitp - halfwidth);
                explogitpup = exp(logitp + halfwidth);
                plo = explogitplo ./ (1 + explogitplo);
                pup = explogitpup ./ (1 + explogitpup);
            end
        end
        
        function y = gampdf(x,a,b)
            %GAMPDF Gamma probability density function.
            %   Y = GAMPDF(X,A,B) returns the gamma probability density function with
            %   shape and scale parameters A and B, respectively, at the values in X.
            %   The size of Y is the common size of the input arguments. A scalar input
            %   functions as a constant matrix of the same size as the other inputs.
            %
            %   Some references refer to the gamma distribution with a single
            %   parameter.  This corresponds to the default of B = 1.
            %
            %   See also GAMCDF, GAMFIT, GAMINV, GAMLIKE, GAMRND, GAMSTAT.
            
            %   References:
            %      [1] Abramowitz, M. and Stegun, I.A. (1964) Handbook of Mathematical
            %          Functions, Dover, New York, section 26.1.
            %      [2] Evans, M., Hastings, N., and Peacock, B. (1993) Statistical
            %          Distributions, 2nd ed., Wiley.
            
            %   Copyright 1993-2004 The MathWorks, Inc.
            %   $Revision: 2.10.2.5 $  $Date: 2004/01/24 09:33:56 $
            
            if nargin < 2
                error('stats:gampdf:TooFewInputs','Requires at least two input arguments');
            elseif nargin < 3
                b = 1;
            end
            
            % Return NaN for out of range parameters.
            a(a <= 0) = NaN;
            b(b <= 0) = NaN;
            
            try
                z = x ./ b;
                
                % Negative data would create complex values, potentially creating
                % spurious NaNi's in other elements of y.  Map them to the far right
                % tail, which will be forced to zero.
                z(z < 0) = Inf;
                
                % Prevent LogOfZero warnings.
                warn = warning('off','MATLAB:log:logOfZero');
                u = (a - 1) .* log(z) - z - gammaln(a);
                warning(warn);
            catch
                if exist('warn','var'), warning(warn); end
                error('stats:gampdf:InputSizeMismatch',...
                    'Non-scalar arguments must match in size.');
            end
            
            % Get the correct limit for z == 0.
            u(z == 0 & a == 1) = 0;
            % These two cases work automatically
            %  u(z == 0 & a < 1) = Inf;
            %  u(z == 0 & a > 1) = -Inf;
            
            % Force a 0 for extreme right tail, instead of getting exp(Inf-Inf)==NaN
            u(z == Inf & isfinite(a)) = -Inf;
            % Force a 0 when a is infinite, instead of getting exp(Inf-Inf)==NaN
            u(z < Inf & a == Inf) = -Inf;
            
            y = exp(u) ./ b;
        end
        
    end
end
