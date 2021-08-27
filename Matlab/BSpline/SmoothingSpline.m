classdef SmoothingSpline < BSpline
    %SmoothingSpline Fit noisy data with a tensioned interpolating spline
    %   3 argument initialization
    %       f = SmoothingSpline(t,x,distribution);
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
        x           % array of values for the dependent axis
        t           % array of values for the independent axis
        distribution % distribution of the noise, must be an instance of a Distribution sublcass.
        
        T           % degree at which tension is applied
        lambda      % smoothing parameter, either pass a numeric value, or the Lambda enumeration. Default is Lambda.optimalIterated
        isConstrained % indicates whether or not lambda was so big that the solution is just a constrained solution
        mu          % mean value of the tension variable
        knot_dof    % knot dofs
        
        covariance  % computed from the given distribution, this is the covariance structure of the observations. It may be a scalar, vector, or matrix.
        
        variableCache % structure storing several cached variables, useful for quick tension spline computation.
        
        sigma       % initial weight (given as normal std dev.)
        
        constraints % constraints = struct('t',[],'D',[]); such that f^(D)(t)=0.
        
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
        CmInv % error in coefficients, MxMxD. This is retrieved from the variable cache.
        nonOutlierIndices
        tensionValue
    end

    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Initialization
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function self = SmoothingSpline(t,x,distribution,varargin)
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
            outlierIndices = [];
            t_knot = [];
            sigma = sqrt(distribution.variance);
            constraints = [];
            
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
                elseif strcmp(varargin{k}, 'sigma')
                    sigma = varargin{k+1};
                elseif strcmp(varargin{k}, 'constraints')
                    constraints = varargin{k+1};
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
                        u_rms = SmoothingSpline.EstimateRMSDerivativeFromSpectrum(t,x,sqrt(distribution.variance),1);
                        n_eff = SmoothingSpline.EffectiveSampleSizeFromUrms(u_rms, t, sqrt(distribution.variance));
                        a_rms = SmoothingSpline.EstimateRMSDerivativeFromSpectrum(t,x,sqrt(distribution.variance),T);
                        lambda = (n_eff-1)/(n_eff*a_rms.^2);
                    case {Lambda.fullTensionExpected, Lambda.optimalIterated}
                        % if you're going to optimize, it's best to start
                        % near the full tension solution, rather than
                        % (potentially) near zero
                        a_rms = SmoothingSpline.EstimateRMSDerivativeFromSpectrum(t,x,sqrt(distribution.variance),T);
                        lambda = 1/a_rms.^2;
                    case  {Lambda.fullTensionIterated, Lambda.optimalExpectedRobust, Lambda.fullTensionExpectedRobust, Lambda.optimalRangedIterated}
                        if length(x)>33
                            x_filtered = SmoothingSpline.RunningFilter(x,11,'median');
                        else
                            x_filtered=x;
                        end
                       a_rms = SmoothingSpline.EstimateRMSDerivativeFromSpectrum(t,x_filtered,sqrt(distribution.variance),T);
                       if lambdaArgument == Lambda.optimalExpectedRobust
                           u_rms = SmoothingSpline.EstimateRMSDerivativeFromSpectrum(t,x_filtered,sqrt(distribution.variance),1);
                           n_eff = SmoothingSpline.EffectiveSampleSizeFromUrms(u_rms, t, sqrt(distribution.variance));
                           lambda = (n_eff-1)/(n_eff*a_rms.^2);
                       else
                           lambda = 1/a_rms.^2;
                       end  
                end
            elseif isscalar(lambdaArgument)
                lambda = lambdaArgument;
            else
                error('Invalid choice for lambda. Lambda must be either a scalar or the enum Lambda.');
            end
                        
            if shouldSetKnotDOFAutomatically == 1
                if isempty(n_eff)
                    u_rms = SmoothingSpline.EstimateRMSDerivativeFromSpectrum(t,x,sqrt(distribution.variance),1);
                    n_eff = SmoothingSpline.EffectiveSampleSizeFromUrms(u_rms, t, sqrt(distribution.variance));
                end
                % conservative estimate
                knot_dof = max(1,floor(0.5*n_eff));
            end
            
            % Compute the spline values at the observation points
            if isempty(t_knot)
                t_knot = InterpolatingSpline.KnotPointsForPoints(t,K,knot_dof);
            end
                        
            if isa(distribution,'NormalDistribution')                
                [m,~,~,isConstrained,cachedVars] = SmoothingSpline.TensionSolution(lambda,mu,t,x,t_knot,K,T,distribution,[],constraints,[]);      
            elseif ~isempty(distribution.w)
                [m,~,~,isConstrained,cachedVars] = SmoothingSpline.IteratedLeastSquaresTensionSolution(lambda,mu,t,x,t_knot,K,T,distribution,constraints,[]);
            else
                error('No weight function given! Unable to proceed.');
            end
                        
            self@BSpline(K,t_knot,m);
            self.lambda = lambda;
            self.isConstrained = isConstrained;
            self.mu = mu;
            self.T = T;
            self.variableCache = cachedVars;
            self.t = t;
            self.x = x;
            
            self.constraints = constraints;
            
            self.knot_dof = knot_dof;
            self.distribution = distribution;
            self.sigma = sigma;
            self.outlierIndices = outlierIndices;
            self.outlierThreshold = self.distribution.locationOfCDFPercentile(1-1/10000/2);
            
            if lambdaArgument == Lambda.optimalIterated
                self.minimizeExpectedMeanSquareError();
            elseif lambdaArgument == Lambda.optimalRangedIterated
                zmin = self.distribution.locationOfCDFPercentile(beta/2);
                zmax = self.distribution.locationOfCDFPercentile(1-beta/2);
                expectedVariance = self.distribution.varianceInRange(zmin,zmax);
                self.minimize( @(spline) spline.expectedMeanSquareErrorInRange(zmin,zmax,expectedVariance) );
            elseif lambdaArgument == Lambda.fullTensionIterated
                self.estimateOutlierDistribution();
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
            if newlambda < 0
                newlambda = 0;
            end
            
            if isempty(self.lambda) 
                self.lambda = newlambda;
            elseif self.lambda ~= newlambda
                self.lambda = newlambda;
                self.tensionParameterDidChange();
            end
        end
        
        function set.outlierIndices(self,outliers)
            self.outlierIndices = outliers;
        end
        
        function nonOutlierIndices = get.nonOutlierIndices(self)
            nonOutlierIndices = setdiff(1:length(self.x),self.outlierIndices);
        end
        
        function self = tensionParameterDidChange(self)
            % Tension parameter changed, so we need to recompute the
            % solution, then then compute the PP coefficients for that
            % solution.
            if isa(self.distribution,'NormalDistribution')
                [self.m,~,~,self.isConstrained,self.variableCache] = SmoothingSpline.TensionSolution(self.lambda,self.mu,self.t,self.x,self.t_knot,self.K,self.T,self.distribution,[],self.constraints,self.variableCache);
            elseif ~isempty(self.distribution.w)
                [self.m,~,~,self.isConstrained,self.variableCache] = SmoothingSpline.IteratedLeastSquaresTensionSolution(self.lambda,self.mu,self.t,self.x,self.t_knot,self.K,self.T,self.distribution,self.constraints,self.variableCache);
            else
                error('No weight function given! Unable to proceed.'); 
            end

            [self.C,self.t_pp,self.B] = BSpline.PPCoefficientsFromSplineCoefficients( self.m, self.t_knot, self.K, self.B );
            
            self.outlierIndices = find(abs(self.epsilon) > self.outlierThreshold);
        end
        
        function set.tensionValue(self,x_T)
            SmoothingSpline.minimizeFunctionOfSpline(self,@(spline) abs(std(spline.uniqueValuesAtHighestDerivative)-x_T));
        end
        
        function x_T = get.tensionValue(self)
           x_T = std(self.uniqueValuesAtHighestDerivative);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Dependent properties
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function X = get.X(self)
            if isempty(self.variableCache.X)
                self.variableCache.X = self.Splines(self.t); 
            end
            X = self.variableCache.X;
        end
        
        function W = get.W(self)
            % The (W)eight matrix used to get the solution. For a normal
            % distribution this is the covariance matrix, but for any other
            % distribution it is not.
            if isempty(self.variableCache.W)
                error('This should not be empty!');
            end
            W = self.variableCache.W;
        end
        
        function XWX = get.XWX(self)
            if isempty(self.variableCache.XWX)
                error('This should not be empty!');
            end
            XWX = self.variableCache.XWX;
        end
        
        function Cm = get.Cm(self)
            if ~isempty(self.variableCache.Cm)
                Cm = self.variableCache.Cm;
            elseif ~isempty(self.variableCache.CmInv)
                Cm = inv(self.variableCache.CmInv);
                self.variableCache.Cm = Cm;
            else
                error('This should not be empty!');
            end
        end
        
        
        function CmInv = get.CmInv(self)
            CmInv = self.variableCache.CmInv;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Smoothing matrix and covariance matrix
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function S = smoothingMatrix(self)
            % The smoothing matrix S takes the observations and maps them
            % onto the estimated true values.
            if size(self.W,1) == length(self.t) && size(self.W,2) == 1
                S = (self.X*(self.CmInv\(self.X.'))).*(self.W.');
            else
                S = (self.X*(self.CmInv\(self.X.')))*self.W;
            end
        end
        
        function S = errorMatrixAtPointsForDerivative(self,t,numDerivs)
            % Returns the covariance matrix for a given derivative at the
            % requested points.
            J = self.Splines(t,numDerivs);            
            if ~isempty(self.distribution.w)
                S = (J*self.Cm*J.')*self.W;
            else
                S = (J*self.Cm*J.')/self.distribution.variance;
            end
        end    
        
        function S = errorMatrixForDerivative(self,numDerivs)
            % Returns the covariance matrix for a given derivative at the
            % points of observation
            S = self.errorMatrixAtPointsForDerivative(self.t,numDerivs);
        end
        
        function SE = standardErrorAtPointsForDerivative(self,t,numDerivs)
            % Returns the standard error for a given derivative at the
            % points requested.
            SE = sqrt(diag(self.errorMatrixAtPointsForDerivative(t,numDerivs)));
        end
        
        function SE = standardErrorForDerivative(self,numDerivs)
            % Returns the standard error for a given derivative at the
            % observation points.
           SE = sqrt(diag(self.errorMatrixForDerivative(numDerivs)));
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % epsilon = observed position minus the fit position
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function epsilon = epsilon(self)
            epsilon = self.x - self.ValueAtPoints(self.t);
        end
        
        function epsilon = epsilonAtIndices(self,indices)
            epsilon = self.x(indices) - self.ValueAtPoints(self.t(indices));
        end
        
        function epsilon = epsilonInRange(self,zmin,zmax)
            epsilon = self.epsilon;
            epsilon = epsilon( epsilon >= zmin & epsilon <= zmax );
        end
        
        function epsilonIQ = epsilonIQ(self)
            % interquartile epsilon
            epsilon = sort(self.epsilon);
            n = length(epsilon);
            epsilonIQ = epsilon(floor(n/4):ceil(3*n/4));
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Sample variance
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function X2 = sampleVariance(self)
            % The sample variance. Same as ||(S-I)*x||^2/N
            % Not normalized by the variance.
            X2 = mean( self.epsilon.^2,1);
        end
        
        function X2 = sampleVarianceInRange(self,zmin,zmax)
            % Sample variance by restricting the computation to epsilons
            % that fall within a certain range.
            epsilon = self.epsilon;
            X2 = mean( epsilon( epsilon >= zmin & epsilon <= zmax ).^2,1);
        end
        
        function X2 = sampleVarianceInPercentileRange(self,pctmin,pctmax)
            zmin = self.distribution.locationOfCDFPercentile(pctmin);
            zmax = self.distribution.locationOfCDFPercentile(pctmax);
            X2 = self.sampleVarianceInRange(zmin,zmax);
        end
        
        function sampleVariance = sampleVarianceIQ(self,alpha)
            % This computes the sample variance from the inner percentage
            % of epsilons, given by alpha. If alpha=1/2 this is the
            % interquartile sample variance.
            if nargin < 2
                alpha = 1/2;
            end
            epsilon = sort(self.epsilon);
            indices = floor(length(epsilon)*alpha/2):ceil(length(epsilon)*(1-alpha/2));
            sampleVariance = mean(epsilon(indices).^2);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Variance of the mean
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function SE2 = varianceOfTheMean(self)
            % The variance of the mean is the square of the standard error
            S = self.smoothingMatrix;
            SE2 = self.distribution.variance*trace(S)/length(S);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Measures of expected error
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
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
        
        function MSE = expectedMeanSquareErrorFromCV(self)
            % Cross-validation (CV) estimate for the mean square error from
            % Green and Silverman, equation 3.5
            epsilon = self.epsilon;
            Sii = diag(self.smoothingMatrix);
            
            MSE = mean( (epsilon./(1-Sii)).^2 );
        end
        
        function MSE = expectedMeanSquareErrorFromGCV(self)
            % Generalized cross-validation (GCV) estimate for the mean
            % square error from Craven and Wahba, 1979 equation 1.9.
            epsilon = self.epsilon;
            S = self.smoothingMatrix;
            
            a = mean(epsilon.^2);
            b = trace(S)/length(S);
            
            MSE = a/(1-b)^2;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Measures of expected error restricted to a subset of points
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [MSE, n] = expectedMeanSquareErrorForPointsAtIndices(self,indices,expectedVariance)
            % This is essentially the same as expectedMeanSquareError
            % above, but restricts the calculation to a subset of points
            % (and let's you set the expected variance).
            epsilon = self.epsilon;
            X2 = mean(epsilon(indices).^2,1);
            
            S = self.smoothingMatrix;
            S = S(indices,indices);
            
            n = length(S);
            MSE = (X2/expectedVariance + 2*trace(S)/n - 1)*expectedVariance;
        end
        
        
        function [MSE, n] = expectedMeanSquareErrorInRange(self,zmin,zmax,expectedVariance)
            epsilon = self.epsilon;
            indices = epsilon >= zmin & epsilon <= zmax;
            if nargin < 4 || isempty(expectedVariance)
                expectedVariance = self.distribution.varianceInRange(zmin,zmax);
            end
            [MSE, n] = self.expectedMeanSquareErrorForPointsAtIndices(indices,expectedVariance);
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
        
        
        function n_eff = effectiveSampleSizeFromSampleVarianceForIndices(self,indices)
            epsilon = self.epsilon;
            n_eff = 1./(1-mean(epsilon(indices).^2)/self.distribution.variance);
        end
        
        function n_eff = effectiveSampleSizeFromVarianceOfTheMeanForIndices(self,indices)
            S = self.smoothingMatrix;
            S = S(indices,indices); 
            n_eff = length(S)/trace(S);
        end
        
        function n_eff = effectiveSampleSizeFromVarianceOfTheMeanInRange(self,zmin,zmax)
            epsilon = self.epsilon;
            n_eff = self.effectiveSampleSizeFromVarianceOfTheMeanForIndices(epsilon >= zmin & epsilon <= zmax);
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
        % Preferred method for achieving full tension
        %
        % 2019/02/26 - The preferred method for now is the iterated
        % interquartile Anderson-Darling test. Other methods are in the
        % ExperimentalSmoothingSpline class
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function estimateOutlierDistribution(self)
            [self.outlierDistribution, self.alpha] = self.setToFullTensionWithIteratedIQAD();
            self.lambdaAtFullTension = self.lambda;
            self.sigmaAtFullTension = sqrt(self.distribution.w(self.epsilon));
        end
        
        function [newOutlierDistribution, newAlpha] = setToFullTensionWithIteratedIQAD(self,knownDistribution)
            % Set to full tension by minimizing the Anderson-Darling (AD)
            % error on the interquartile (IQ) range of epsilons. At each
            % iteration the outlier distribution is estimated and used to
            % refine the tension.
            if nargin < 2
                knownDistribution = self.distribution;
            end
            
            self.minimize(@(spline) knownDistribution.andersonDarlingInterquartileError(spline.epsilon));
            lastAlpha = 0.0;
            totalIterations = 0;
            [newOutlierDistribution, newAlpha] = SmoothingSpline.estimateOutlierDistributionFromKnownNoise(self.epsilon,knownDistribution);
            % The estimateOutlierDistributionFromKnownNoise function always
            % works from the same list of 100 alphas, so we can check with
            % fairly high tolerance whether or not its the same.
            while (abs(lastAlpha-newAlpha) > 1e-6 && totalIterations < 10)
                if newAlpha > 0 && ~isempty(newOutlierDistribution)
                    addedDistribution = AddedDistribution(newAlpha,newOutlierDistribution,knownDistribution);
                else
                    addedDistribution = knownDistribution;
                end
                self.minimize(@(spline) addedDistribution.andersonDarlingInterquartileError(spline.epsilon));
                lastAlpha = newAlpha;
                [newOutlierDistribution, newAlpha] = SmoothingSpline.estimateOutlierDistributionFromKnownNoise(self.epsilon,knownDistribution);
                totalIterations = totalIterations + 1;
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Minimization
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [lambda,fval] = minimize(self,penaltyFunction)
            % the penalty function *must* take a tension spline object and
            % return a scalar. The function will be minimized by varying
            % lambda.
            [~,lambdaBelow,lambdaAbove] = SmoothingSpline.minimizeFunctionOfSplineWithGridSearch(self,penaltyFunction);
            lambda = SmoothingSpline.minimizeFunctionOfSplineBounded(self,penaltyFunction,lambdaBelow,lambdaAbove);
            fval = penaltyFunction(self);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Minimization convenience methods
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [lambda,mse] = minimizeMeanSquareError(self,t_true,x_true)
            % minimize against the true values
            mse = @(aSmoothingSpline) mean((aSmoothingSpline.ValueAtPoints(t_true)-x_true).^2)/(aSmoothingSpline.distribution.variance);
            [lambda,mse] = self.minimize(mse);
        end
        
        function [lambda,mse] = minimizeExpectedMeanSquareError(self)
            % minimize against the *expected* mean square error for the
            % given distribution
            [lambda,mse] = self.minimize(@(aSmoothingSpline) aSmoothingSpline.expectedMeanSquareError);
        end
        
        function [lambda,mse] = minimizeExpectedMeanSquareErrorInPercentileRange(self,pctmin,pctmax)
            zmin = self.distribution.locationOfCDFPercentile(pctmin);
            zmax = self.distribution.locationOfCDFPercentile(pctmax); 
            expectedVariance = self.distribution.varianceInRange(zmin,zmax);
            [lambda,mse] = self.minimize( @(aSmoothingSpline) aSmoothingSpline.expectedMeanSquareErrorInRange(zmin,zmax,expectedVariance) );
        end
        
        function [lambda,mse] = minimizeExpectedMeanSquareErrorInNoiseRange(self,minimizationPDFRatio)
            if nargin < 2
                minimizationPDFRatio = 1;
            end
            if ~isempty(self.alpha) && self.alpha > 0
                [zoutlier,expectedVariance] = SmoothingSpline.locationOfNoiseToOutlierPDFRatio(self.alpha,self.outlierDistribution,self.distribution,minimizationPDFRatio);
                [lambda,mse] = self.minimize( @(spline) spline.expectedMeanSquareErrorInRange(-zoutlier,zoutlier,expectedVariance));
            else
                [lambda,mse] = self.minimizeExpectedMeanSquareError();
            end
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
        %
        % This algorithm also *stops* looking at higher lambdas once the
        % problem becomes constrained.
        function [lambda,lambdaBelow,lambdaAbove] = minimizeFunctionOfSplineWithGridSearch(aSmoothingSpline,functionOfSpline,dlog10lambda,n)
            lambdas = aSmoothingSpline.lambda;
            err = functionOfSpline(aSmoothingSpline);
            isConstrained = aSmoothingSpline.isConstrained;
            nLambdas = 1;
            index = 1;
            
            if nargin < 4
                dlog10lambda = 0.5; % step size in factors of 10
                n = 3; % number of steps to make each direction
            end
            
            while (index <= 2 || (index >= nLambdas-1 && isConstrained(index) ~=1) )   
                if index <= 2
                    % expand search below
                    lambda0 = lambdas(1);
                    lambdas = cat(1,10.^linspace(log10(lambda0)-n*dlog10lambda,log10(lambda0)-dlog10lambda,n)',lambdas);
                    err = cat(1,zeros(n,1),err);
                    isConstrained = cat(1,zeros(n,1),isConstrained);
                    for iLambda = 1:n
                        aSmoothingSpline.lambda = lambdas(iLambda);
                        err(iLambda) =  functionOfSpline(aSmoothingSpline);
                        isConstrained(iLambda) = aSmoothingSpline.isConstrained;
                    end
                    
                    nLambdas = length(lambdas);
                    index = index+n;
                end
                
                if (index >= nLambdas-1 && isConstrained(index) ~=1)
                    % expand search above
                    lambda0 = lambdas(end);
                    lambdas = cat(1,lambdas,10.^linspace(log10(lambda0)+dlog10lambda,log10(lambda0)+n*dlog10lambda,n)');
                    err = cat(1,err,zeros(n,1));
                    isConstrained = cat(1,isConstrained,zeros(n,1));
                    for iLambda = (nLambdas+1):(nLambdas+n)
                        aSmoothingSpline.lambda = lambdas(iLambda);
                        err(iLambda) =  functionOfSpline(aSmoothingSpline);
                        isConstrained(iLambda) = aSmoothingSpline.isConstrained;
                    end
                    
                    nLambdas = length(lambdas);
                end
                
                [minErr,index] = min(err);
                if index <= 2
                    % if we reach a plateau at the bottom, we need to stop
                    index = find(err<=minErr,1,'last');
                end
            end
            lambda = lambdas(index);
            aSmoothingSpline.lambda = lambda; 
            lambdaBelow = lambdas(index-1);
            if length(lambdas)>index
                lambdaAbove = lambdas(index+1);
            else
                lambdaAbove = lambdas(index);
            end
        end
        
        function lambda = minimizeFunctionOfSplineBounded(aSmoothingSpline,functionOfSpline,x1,x2)
            epsilon = 1e-15;
            errorFunction = @(log10lambdaPlusEpsilon) SmoothingSpline.FunctionOfSplineWrapper(aSmoothingSpline,log10lambdaPlusEpsilon,functionOfSpline);
            optimalLog10lambdaPlusEpsilon = fminbnd( errorFunction, log10(x1+epsilon), log10(x2+epsilon), optimset('TolX', 0.01, 'TolFun', 0.001) );
            lambda = 10^optimalLog10lambdaPlusEpsilon - epsilon;
            aSmoothingSpline.lambda = lambda;
        end
        
        function lambda = minimizeFunctionOfSpline(aSmoothingSpline,functionOfSpline)
            epsilon = 1e-15;
            errorFunction = @(log10lambdaPlusEpsilon) SmoothingSpline.FunctionOfSplineWrapper(aSmoothingSpline,log10lambdaPlusEpsilon,functionOfSpline);
            optimalLog10lambdaPlusEpsilon = fminsearch( errorFunction, log10(aSmoothingSpline.lambda+epsilon), optimset('TolX', 0.01, 'TolFun', 0.001) );
            lambda = 10^optimalLog10lambdaPlusEpsilon - epsilon;
            aSmoothingSpline.lambda = lambda;
        end
        
        function error = FunctionOfSplineWrapper(aSmoothingSpline, log10lambdaPlusEpsilon, functionOfSpline)
            epsilon = 1e-15;
            aSmoothingSpline.lambda = 10^log10lambdaPlusEpsilon-epsilon;              
            error = functionOfSpline(aSmoothingSpline);
        end
        
        % Same thing as above, but for a bunch of splines given as a cell
        % array
        function lambda = minimizeFunctionOfSplines(tensionSplines,functionOfSplines)
            epsilon = 1e-15;
            errorFunction = @(log10lambdaPlusEpsilon) SmoothingSpline.functionOfSplinesWrapper(tensionSplines,log10lambdaPlusEpsilon,functionOfSplines);
            optimalLog10lambdaPlusEpsilon = fminsearch( errorFunction, log10(tensionSplines{1}.lambda+epsilon), optimset('TolX', 0.01, 'TolFun', 0.001) );
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

        function cachedVars = PrecomputeTensionSolutionMatrices(t,x,t_knot,K,T,distribution,W,constraints,cachedVars)
            cachedVars = ConstrainedSpline.PrecomputeSolutionMatrices(t,x,K,t_knot,distribution,W,constraints,cachedVars);
            
            if ~isfield(cachedVars,'m_constrained') || ~isfield(cachedVars,'Cm_constrained')
                cachedVars.m_constrained = [];
                cachedVars.Cm_constrained = [];
            end
            
            N = length(x);
            if ~isfield(cachedVars,'V') || isempty(cachedVars.V)
                % This is the value of the tensioned variable on the
                % quadrature grid.
                Q = 10*N; % number of points on the quadrature grid
                tq = linspace(t(1),t(end),Q)';
                B = BSpline.Spline( tq, t_knot, K, T );
                cachedVars.V = squeeze(B(:,:,T+1)); % QxM
            end
            
            if ~isfield(cachedVars,'VV') || isempty(cachedVars.VV)
                V = cachedVars.V;
                cachedVars.VV = V'*V;
            end
            
 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Repeat above, but now with shape constraints
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            S = cachedVars.SC;
            if ~isempty(S)    
                if ~isfield(cachedVars,'VS') || isempty(cachedVars.VS)
                    cachedVars.VS = cachedVars.V*S;
                end        
                
                if ~isfield(cachedVars,'SVVS') || isempty(cachedVars.SVVS)
                    VS = cachedVars.VS;
                    cachedVars.SVVS = VS'*VS;
                end
            end
        end
        
        function [m,Cm,CmInv,isConstrained,cachedVars] = TensionSolution(lambda,mu,t,x,t_knot,K,T,distribution,W,constraints,cachedVars)
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
            
            cachedVars = SmoothingSpline.PrecomputeTensionSolutionMatrices(t,x,t_knot,K,T,distribution,W,constraints,cachedVars);
            
            F = cachedVars.F;
            XWX = cachedVars.XWX;
            XWx = cachedVars.XWx;
            VV = cachedVars.VV;
            
            N = length(x); % number of points
            M = size(XWX,2); % number of splines
            Q = size(cachedVars.V,1); % number of quadrature points
            
            % set up inverse matrices
            E_x = XWX + (lambda*N/Q)*(VV);
            F_x = XWx;
            
            % add the mean tension value
            if mu ~= 0.0
                V = cachedVars.V;
                F_x = F_x + (lambda*N/Q)*mu*transpose(sum( V,1));
            end
            
            % Expand with local constraints
            NC = size(F,1);
            if NC > 0
                E_x = cat(1,E_x,F); % (M+NC)xM
                E_x = cat(2,E_x,cat(1,F',zeros(NC)));
                F_x = cat(1,F_x,zeros(NC,1));
            end
            
            lastwarn('');
            m = E_x\F_x;
            [~,msgID] = lastwarn;
            
            if strcmp(msgID,'MATLAB:nearlySingularMatrix')
                % Check if the tension is sufficiently high that we just need
                % to use the constrained solution. This scenario only makes
                % sense because we know we've used the interpolation points as
                % knot points... otherwise rcond->0 could mean something else.
                if isempty(cachedVars.m_constrained) || isempty(cachedVars.Cm_constrained)
                    
                    t_knot_constrained = cat(1,min(t)*ones(T,1),max(t)*ones(T,1));
                    % T=2 indicates tension on acceleration, which we want
                    % to be zero, so we would want K=2
                    cspline = ConstrainedSpline(t,x,T,t_knot_constrained,distribution,constraints);
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
                
                if NC > 0
                    M = size(cachedVars.X,2);
                    m = m(1:M);                    
                end
                
                S = cachedVars.SC;
                if ~isempty(S)
                    m0 = m;
                    xi0 = S\m0;
                    if any(xi0<0)
                        E_x = cachedVars.SXWXS + (lambda*N/Q)*(cachedVars.SVVS); % MxM
                        F_x = cachedVars.SXWx;
                        
                        lb = zeros(M,1); lb(1) = -inf;
                        ub = inf*ones(M,1);
                        
                        % These are the local constraints, exactly as above
                        Aeq = cachedVars.FS;
                        if isempty(Aeq)
                            beq = [];
                        else
                            beq = zeros(NC,1);
                        end
                        
                        % E_x should be symmetric, although sometimes it's not
                        % exactly.
                        H = (E_x+E_x')*0.5;
                        
                        if NC == 0
                            options = optimoptions('quadprog','Display','off','Algorithm','trust-region-reflective');
                            x = quadprog(2*H,-2*F_x,[],[],Aeq,beq,lb,ub,xi0,options);
                        else
                            options = optimoptions('quadprog','Display','off','Algorithm','interior-point-convex');
                            x = quadprog(2*H,-2*F_x,[],[],Aeq,beq,lb,ub,[],options);
                        end
                        m = S*x;
                    end
                end
                
                Cm = [];
                if NC > 0
                    CmInv = XWX + (lambda*N/Q)*(VV);
                else
                    CmInv = E_x;
                end
                
            end
            
            cachedVars.Cm = Cm;
            cachedVars.CmInv = CmInv;
        end
        
        % e
        function [m,Cm,CmInv,isConstrained,cachedVars] = IteratedLeastSquaresTensionSolution(lambda,mu,t,x,t_knot,K,T,distribution,constraints,cachedVars)
            % Out first call is basically just a normal fit to a tension
            % spline. If W is not set, it will be set to either 1/w0^2 or
            % the correlated version
            cachedVars.W = []; cachedVars.XWX = []; cachedVars.XWx = [];
            [m,Cm,CmInv,isConstrained,cachedVars] = SmoothingSpline.TensionSolution(lambda,mu,t,x,t_knot,K,T,distribution,[],constraints,cachedVars);
            
            X = cachedVars.X;
            sigma2_previous = (distribution.sigma0)^2;
            rel_error = 1.0;
            repeats = 1;
            while (rel_error > 0.01)
                sigma_w2 = distribution.w(X*m - x);
                
                if exist('cachedVars.rho_t','var')
                    Sigma2 = (sqrt(sigma_w2) * sqrt(sigma_w2).') .* cachedVars.rho_t;
                    W = inv(Sigma2);
                else
                    W = 1./sigma_w2;
                end
                
                % hose any cached variable associated with W...
                cachedVars.W = []; cachedVars.XWX = []; cachedVars.XWx = [];
                % ...and recompute the solution with this new weighting
                [m,Cm,CmInv,isConstrained,cachedVars] = SmoothingSpline.TensionSolution(lambda,mu,t,x,t_knot,K,T,distribution,W,constraints,cachedVars);
                
                rel_error = max( abs((sigma_w2-sigma2_previous)./sigma_w2) );
                sigma2_previous=sigma_w2;
                repeats = repeats+1;
                
                if (repeats == 250)
                    disp('Failed to converge after 250 iterations.');
                    break;
                end
            end
            
        end
         
        function [m,Cm] = CoefficientsForConstrainedSolution(t,x,K,distribution,X)
            % The constrained solution is the solution where (d/dt)^(T) =
            % 0, which is the same as a K=T spline.
            t_knot = cat(1,min(t)*ones(K,1),max(t)*ones(K,1));
            cspline = ConstrainedSpline(t,x,K,t_knot,distribution,[]);
            m = X\cspline(t);
            
            % if m_x = inv(X)*Y*m_y then,
            % Cm_x = G*Cm_y*G^T if G = inv(X)*Y
            Y = cspline.X;
            G = X\Y;
            Cm = G*(CmyInv\(G.'));          
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Supporting function
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [sigma, mu] = StandardDeviationAndMeanOfInterquartileRange(x)
            % am I sure this is right? This looks wrong!
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
        % Searching the distribution pdfs for specific points
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [z,expectedVariance] = locationOfNoiseToOutlierPDFRatio(alpha,outlierDistribution,noiseDistribution,pdfRatio)
            % Setting the pdfRatio to 1 means that this will return the
            % point where the odds of being an outlier are equal with the
            % odds of being characterized noise. 1 seems to be a good
            % value... but a bit higher is good too, like 3 or so. This
            % value is *less* than the variance crossover found with the
            % above distribution.
            f = @(z) abs((1-alpha)*noiseDistribution.pdf(z)./(alpha*outlierDistribution.pdf(z))-pdfRatio);
            z = abs(fminsearch(f,sqrt(noiseDistribution.variance)));
            expectedVariance = AddedDistribution(alpha,outlierDistribution,noiseDistribution).varianceInRange(-z,z);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Estimating the outleir distribution
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [outlierDistribution, alpha] = estimateOutlierDistributionFromKnownNoise(epsilon,noiseDistribution)
            % Given epsilon at full tension, and some characterized noise
            % distribution, estimate the outlier distribution.
            %
            % This method generates all possible fractions "added
            % distributions" that enough variance to match the observed
            % value, and then find the fraction that minimizes the
            % Anderson-Darling test.
            
            % Let's find a reasonable set of z_crossover points.
            s2_total = mean(epsilon.^2);
            s2_noise = noiseDistribution.variance;
            
            % What's a reasonable cutoff here? For now 2.0 is just a magic
            % number...
            if s2_total/s2_noise < 2.0
                outlierDistribution = [];
                alpha = 0;
                return;
            end
            
            alpha_outlier = 10.^(linspace(log10(0.01),log10(0.5),100))';
            sigma2_outlier = (s2_total-(1-alpha_outlier)*s2_noise)./(3*alpha_outlier);
            var_total = zeros(size(alpha_outlier));
            ks_error = zeros(size(alpha_outlier));
            
            for iAlpha = 1:length(alpha_outlier)
                newAddedDistribution = AddedDistribution(alpha_outlier(iAlpha),StudentTDistribution(sqrt(sigma2_outlier(iAlpha)),3),noiseDistribution);
                var_total(iAlpha) = newAddedDistribution.variance;
                ks_error(iAlpha) = newAddedDistribution.andersonDarlingError(epsilon);
            end
            
            [minError,minIndex] = min(ks_error);
            if (noiseDistribution.andersonDarlingError(epsilon) < minError)
                outlierDistribution = [];
                alpha = 0;
            else
                alpha = alpha_outlier(minIndex);
                sigma2o = sigma2_outlier(minIndex);
                outlierDistribution = StudentTDistribution(sqrt(sigma2o),3);
            end
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Methods for estimating lambda
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
        function n_eff = EffectiveSampleSizeFromUrms(u_rms, t, sigma)
            % These are the coefficients of the empirical best fits for
            % slopes [-2,-3,-4] to the model n_eff = exp(b)*gamma^m
            m = [0.6935; 0.7147; 0.7016];
            b = [2.2594; 2.6123; 2.7663];
            
            dt = median(abs(diff(t)));
            
            gamma = sigma/( u_rms*dt );
            n_eff = max(1,exp(b(2))*gamma.^m(2));
        end
                        
        function [a_rms, a_std, a_mean] = EstimateRMSDerivativeFromSpectrum( t, x, sigma, D, shouldPlotSpectra)
            % EstimateRMSDerivativeFromSpectrum Given some signal (t,x)
            % contaminated by noise sigma, this uses the spectrum to
            % estimate u_rms.
            %
            % D = 1 is velocity, D=2 is acceleration
            
            xin = x;
            tin = t;
            
            if SmoothingSpline.IsEvenlySampled(t) ~= 1
                %    fprintf('interpolating...\n');
                dt = median(diff(t));
                N = ceil((t(end)-t(1))/dt);
                t2 = dt*((0:(N-1))') + t(1);
                x = interp1(t,x,t2);
                t = t2;
            end
            
            [p,~,mu]=polyfit(t,x,D);
            a_mean = factorial(D)*p(1)/mu(2)^D;
            
            % now remove the trend
            x = x-polyval(p,t,[],mu);
                  
            if 1 == 1
                
                dt = t(2) - t(1);
                % temp hack?begin
                t(end) = [];
                x = diff(x)/dt;
                D = D-1;
                % temp hack?end
                
                T = t(end)-t(1);
                N = length(t);
                
                df = 1/T;
                f = ([0:ceil(N/2)-1 -floor(N/2):-1]*df)';
                
                SpectralD = (2*pi*f).^(2*D);
                SpectralD = (2*(1-cos(dt*2*pi*f))/(dt^2)).^D;
                
                ubar = fft(x);
                s_signal = (ubar.*conj(ubar)) .* SpectralD * (dt/N);
                
                SpectralD = (2*(1-cos(dt*2*pi*f))/(dt^2)).^(D+1);
            else
                [DiffMatrix,t_u] = SmoothingSpline.FiniteDifferenceMatrixNoBoundary(D,t,1);
                
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
            s_noise = sigma*sigma*dt*SpectralD;
            
            % The factor of 10 is consistent with 80% confidence.
            % 95% confidence (actually, 97.5% ?) is 39.5
            alpha = 0.10;
            cutoff = 2/SmoothingSpline.chi2inv(alpha/2,2);
            
            % There are two ways to think of this. Either you look for the
            % 95% range of the signal, or the 95% range of the expected
            % noise.
            alpha = 0.99999;
            dof = 2;
            cutoff = 1*SmoothingSpline.chi2inv(alpha,dof)/dof;
            
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
                
                figure
                plot(tin,xin), hold on, plot(tin,polyval(p,tin,[],mu))
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
        
        
        function [rmse,norm] = MeanSquareErrorAtAllOrders(aSmoothingSpline, t_true, x_true)
            if length(unique(diff(t_true))) > 1
                error('This only works for evenly spaced data at the moment.');
            end
            
            rmse = zeros(aSmoothingSpline.K,1);
            norm = zeros(aSmoothingSpline.K,1);
            dt = t_true(2)-t_true(1);
            for D = 0:(aSmoothingSpline.K-1)
                
                % differentiate D times
                % We remove the mean from *position*
                signal_true = x_true;
                for i=1:D
                    signal_true = diff(signal_true)/dt;
                end
                
                if 1 == 0
                     % points of accuracy move dt/2 further inside each D
                    t_signal = t_true(1:(length(t_true)-D)) + D*dt/2;
                    
                    signal_obs = aSmoothingSpline.ValueAtPoints(t_signal,D);
                else
                    signal_obs = aSmoothingSpline.ValueAtPoints(t_true);
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
        
        function Q = QScore(aSmoothingSpline, t_true, x_true)
            [rmse,norm] = SmoothingSpline.MeanSquareErrorAtAllOrders(aSmoothingSpline, t_true, x_true);
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
                
                c = SmoothingSpline.weights( z(i), x(range), numDerivs );
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
        % Smoothing
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Computes a running average across ensembles and reports standard error
        % size(vec) = [n m] where n is time, and m is the ensemble
        
        function [y,yerr] = RunningFilter( x, smoothness, filter )
            %	Returns a running average of vec. Tries to handle end points well.
            
            n = size(x,1);
            m = size(x,2);
            
            y = NaN(n,1);
            yerr = NaN(n,1);
            
            smoothnessHalfLength = ceil((smoothness - 1)/2);
            
            for index=1:n
                restrictionDistance = min([smoothnessHalfLength; index-1; n-index]);
                indices = (index-restrictionDistance):(index+restrictionDistance);
                
                if ( ~isempty(indices) )
                    vals = reshape(x(indices,:),[length(indices)*m 1]);
                    if strcmp(filter,'mean')
                        y(index) = mean(vals);
                    elseif strcmp(filter,'median')
                        y(index) = median(vals);
                    end
                    yerr(index) = std(vals)/sqrt(length(indices)*m);
                end
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
            x = SmoothingSpline.gaminv(p,v/2,2);
            
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
                
                F = SmoothingSpline.gamcdf(q,a,1);
                f = max(SmoothingSpline.gampdf(q,a,1), realmin(class(p)));
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
