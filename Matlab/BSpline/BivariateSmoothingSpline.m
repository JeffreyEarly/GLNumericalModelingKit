classdef BivariateSmoothingSpline < handle
    % BivariateSmoothingSpline Fit 2D (two-dependent, one-independent) data
    % isotropically.
    %
    %   4 argument initialization
    %       f = BivariateSmoothingSpline(t,x,y,distribution);
    %   where
    %       t               time (seconds or datetime)
    %       x               latitude
    %       y               longitude
    %       distribution    distribution of errors
    %       f               spline interpolant
    %
    %   BivariateSmoothingSpline takes a number of optional input argument pairs.
    
    properties
        x
        y
        t
        distribution
        
        K           % spline order (degree+1)
        T           % degree at which tension is applied
        lambda      % tension parameter
        
        x_bar       % x mean (low passed, at observations times)
        y_bar       % y mean (low passed, at observations times)
        x_prime     % x variation (observations, minus low pass)
        y_prime     % y variation (observations, minus low pass)
        
        shouldUseRobustFit = 0;
        % These are not used directly in minimization, only the
        % 'distribution' property in the SmoothingSpline superclass is used
        % directly.
        noiseDistribution
        outlierDistribution
        alpha
        noiseDistanceDistribution
        outlierDistanceDistribution
        sigma
        
        % These have no consequence to the fit, but can be useful for
        % diagnostics.
        outlierIndices = [] 
        outlierThreshold % set to a distance with < 1/10000 odds
        
        lambdaAtFullTension % what was the value of lambda at full tension
        sigmaAtFullTension % what was the set of 'sigmas' produced from the full tension solution
        
        % lookup table for tension
        tensionLambdaTable = [] % size(n,3) [x_T lambda_x lambda_y];
        
        % splines
        spline_mean_x
        spline_mean_y
        spline_x
        spline_y
    end
    
    properties (Dependent)
        tensionValue
    end
    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Initialization
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function self = BivariateSmoothingSpline(t,x,y,distribution,varargin)
            N = length(t);
            t = reshape(t,[],1);
            x = reshape(x,[],1);
            y = reshape(y,[],1);
            
            if length(x) ~= N || length(y) ~= N
                error('x, y, and t must have the same length.');
            end
            
            if ~isa(distribution,'Distribution')
                error('The distribution must be a Distribution subclass.')
            end
            
            nargs = length(varargin);
            if mod(nargs,2) ~= 0
                error('Arguments must be given as name/value pairs');
            end
            
            self.K = 4; % default spline order (degree+1)
            self.T = self.K-1; % default tension *degree* (order-1)
            
            for k = 1:2:length(varargin)
                if strcmp(varargin{k}, 'K')
                    self.K = varargin{k+1};
                elseif strcmp(varargin{k}, 'S')
                    self.K = varargin{k+1}+1;
                elseif strcmp(varargin{k}, 'T')
                    self.T = varargin{k+1};
                elseif strcmp(varargin{k}, 'shouldUseRobustFit')
                    self.shouldUseRobustFit = varargin{k+1};
                end
            end
            
            self.x = x;
            self.y = y;
            self.t = t;
            self.distribution = distribution;
            self.noiseDistribution = distribution;
            self.noiseDistanceDistribution = TwoDimDistanceDistribution(self.noiseDistribution);
            self.outlierThreshold = self.noiseDistanceDistribution.locationOfCDFPercentile(1-1/10000/2);
            
            t_knot = cat(1,min(t)*ones(self.K+1,1),max(t)*ones(self.K+1,1));
            self.spline_mean_x = ConstrainedSpline(self.t,self.x,self.K+1,t_knot,self.distribution,[]);
            self.spline_mean_y = ConstrainedSpline(self.t,self.y,self.K+1,t_knot,self.distribution,[]);
            
            self.x_bar = self.spline_mean_x(self.t);
            self.y_bar = self.spline_mean_y(self.t);
            
            self.x_prime = self.x - self.x_bar;
            self.y_prime = self.y - self.y_bar;
             
            % initialize splines so that they have the same lambda
            if self.shouldUseRobustFit == 1
                % We don't let the RobustSmoothingSpline do any extra work
%                 nu = 3.0;
%                 sigma = sqrt(self.noiseDistribution.variance*1000*(nu-2)/nu);
%                 addedDistribution = AddedDistribution(1/100,StudentTDistribution(sigma,nu),self.noiseDistribution);
                self.spline_x = SmoothingSpline(self.t,self.x_prime,self.noiseDistribution,'K',self.K,'T',self.T,'lambda',Lambda.fullTensionIterated);
                self.spline_y = SmoothingSpline(self.t,self.y_prime,self.noiseDistribution,'K',self.K,'T',self.T,'lambda',self.spline_x.lambda);
                self.lambda = self.spline_x.lambda;
%                 
%                 self.spline_x.sigma = self.noiseDistribution.sigma;
%                 self.spline_y.sigma = self.noiseDistribution.sigma;
                
                % Step 2---minimize under current conditions
%                 [mse1_lambda,mse1] = self.minimizeExpectedMeanSquareErrorInNoiseRange();
                [mse1_lambda,mse1] = self.minimizeExpectedMeanSquareErrorInPercentileRange(1-1/100);
                
                return
                
                % Step 1---estimate the outlier distribution
                self.estimateOutlierDistribution();
                self.lambdaAtFullTension = self.lambda;
                self.sigmaAtFullTension = sqrt(self.spline_x.distribution.w(self.epsilon_d/sqrt(2)));

                self.spline_x.sigma = self.noiseDistribution.sigma;
                self.spline_y.sigma = self.noiseDistribution.sigma;
                
                % Step 2---minimize under current conditions
%                 [mse1_lambda,mse1] = self.minimizeExpectedMeanSquareErrorInNoiseRange();
                [mse1_lambda,mse1] = self.minimizeExpectedMeanSquareErrorInPercentileRange(1-1/100);
                return
                % Step 3---apply outlier method
                self.spline_x.sigma = self.sigmaAtFullTension;
                self.spline_y.sigma = self.sigmaAtFullTension;
                
                % Step 4---re-compute minimum expected mean square error
%                 [~,mse2] = self.minimizeExpectedMeanSquareErrorInNoiseRange();
                [~,mse2] = self.minimizeExpectedMeanSquareErrorInPercentileRange(1-1/100);
                
                % Step 5---set to global minimum
                if mse1 < mse2
                    self.spline_x.sigma = sqrt(self.spline_x.distribution.variance);
                    self.spline_y.sigma = sqrt(self.spline_x.distribution.variance);
                    self.lambda = mse1_lambda;
                end
            else
                self.spline_x = SmoothingSpline(self.t,self.x_prime,self.distribution,'K',self.K,'T',self.T);
                self.spline_y = SmoothingSpline(self.t,self.y_prime,self.distribution,'K',self.K,'T',self.T,'lambda',self.spline_x.lambda);
                self.lambda = self.spline_x.lambda;
                
                self.minimizeExpectedMeanSquareErrorInNoiseRange();
            end
            
            
            
        end
        
%         function [Sx, Sy] = smoothingMatrices(self)
%             Sbar = self.spline_mean_x.smo
%         end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Evaluation
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [x,y] = xyAtTime(self,time)
            x = self.spline_x(time) + self.spline_mean_x(time);
            y = self.spline_y(time) + self.spline_mean_y(time);
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
        
        function self = tensionParameterDidChange(self)
            self.spline_x.lambda = self.lambda;
            self.spline_y.lambda = self.lambda;
            self.outlierIndices = find(self.epsilon_d > self.outlierThreshold);
        end
                
                
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % epsilon = observed position minus the fit position
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function epsilon = epsilon(self)
            epsilon = cat(2,self.spline_x.epsilon,self.spline_y.epsilon);
        end
        
        function epsilon_d = epsilon_d(self)
            epsilon_d = sqrt(sum(self.epsilon.^2,2));
        end
        
        function a = isConstrained(self)
           a = self.spline_x.isConstrained && self.spline_y.isConstrained;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % smoothing matrix
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function Sx = smoothingMatrixX(self)
            Sbar = self.spline_mean_x.smoothingMatrix;
            Sprime = self.spline_x.smoothingMatrix;
            Sx = Sbar + Sprime - Sprime*Sbar;
        end
        
        function Sy = smoothingMatrixY(self)
            Sbar = self.spline_mean_y.smoothingMatrix;
            Sprime = self.spline_y.smoothingMatrix;
            Sy = Sbar + Sprime - Sprime*Sbar;
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Measures of expected error
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function MSE = expectedMeanSquareError(self)
            % This is the *expected* mean-square error normalized by the
            % variance, simply added together across dimensions. Does *not*
            % include the mean.
            
            MSE = self.spline_x.expectedMeanSquareError + self.spline_y.expectedMeanSquareError;
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
            n=sum(indices>0);
            MSE = self.spline_x.expectedMeanSquareErrorForPointsAtIndices(indices,expectedVariance) + self.spline_x.expectedMeanSquareErrorForPointsAtIndices(indices,expectedVariance);
        end
        
        function [MSE, n] = expectedMeanSquareErrorInDistanceRange(self,zmin,zmax,expectedVariance)
            epsilon = self.epsilon_d;
            indices = epsilon >= zmin & epsilon <= zmax;
            if nargin < 4 || isempty(expectedVariance)
                expectedVariance = self.distribution.varianceInRange(zmin,zmax);
            end
            [MSE, n] = self.expectedMeanSquareErrorForPointsAtIndices(indices,expectedVariance);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % 2D methods for achieving full tension
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function setToFullTensionWithIteratedIQAD(self)
            % Set to full tension by minimizing the Anderson-Darling (AD)
            % error on the interquartile (IQ) range of epsilons. At each
            % iteration the outlier distribution is estimated and used to
            % refine the tension.
            if self.alpha > 0
                initialDistribution = AddedDistribution(self.alpha,self.outlierDistribution,self.distribution);
            else
                initialDistribution = self.distribution;
            end
            self.minimize(@(spline) initialDistribution.andersonDarlingInterquartileError( reshape(spline.epsilon,[],1) ));
            
            lastAlpha = 0.0;
            totalIterations = 0;
            [newOutlierDistribution, newAlpha] = RobustSmoothingSpline.estimateOutlierDistributionFromKnownNoise(reshape(self.epsilon,[],1),self.distribution);
            while (abs(lastAlpha-newAlpha) > 0.01 && totalIterations < 10)
                if newAlpha > 0 && ~isempty(newOutlierDistribution)
                    addedDistribution = AddedDistribution(newAlpha,newOutlierDistribution,self.distribution);
                else
                    addedDistribution = self.distribution;
                end
                self.minimize(@(spline) addedDistribution.andersonDarlingInterquartileError( reshape(spline.epsilon,[],1) ));
                lastAlpha = newAlpha;
                [newOutlierDistribution, newAlpha] = RobustSmoothingSpline.estimateOutlierDistributionFromKnownNoise(reshape(self.epsilon,[],1),self.distribution);
                if newAlpha >= 0.5
                    fprintf('Alpha reached 0.5!\n'); break;
                end
                totalIterations = totalIterations + 1;
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Methods for outlier identification
        %
        % We include several strategies. All strategies require estimating
        % outlier distribution to determine the outliers, and the first
        % three strategies use some sort of ranged minimization.
        %
        %   1. Prime the IRLS with a new sigma, tuned to outliers,
        %   2. Change the IRLS weight function for outliers,
        %   3. Remove knot support for outliers,
        %   4. Weight the values of epsilon according to their outlier
        %   probability, and used a weighted MSE.
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Note: I could universally have a 'full tension' function that
        % identifies outliers as the pdf ratio for ALL these functions.
        
        % Minimization should be separated out in this same way.
        
        function estimateOutlierDistribution(self)
            self.setToFullTensionWithIteratedIQAD();
            
            [self.outlierDistribution,self.alpha] = RobustSmoothingSpline.estimateOutlierDistributionFromKnownNoise(reshape(self.epsilon,[],1),self.noiseDistribution);
            if self.alpha>0
                self.outlierDistanceDistribution = TwoDimDistanceDistribution(self.outlierDistribution);
            end
        end
        
        function setSigmaFromFullTensionSolution(self)
            self.setToFullTensionWithIteratedIQAD();
            
            [self.outlierDistribution,self.alpha] = RobustSmoothingSpline.estimateOutlierDistributionFromKnownNoise(reshape(self.epsilon,[],1),self.noiseDistribution);
            self.noiseDistanceDistribution = TwoDimDistanceDistribution(self.noiseDistribution);
            self.outlierDistanceDistribution = TwoDimDistanceDistribution(self.outlierDistribution);
            
            self.sigma = self.noiseDistribution.w(self.epsilon_d/sqrt(2));
            self.spline_x.sigma = self.sigma;
            self.spline_y.sigma = self.sigma;
        end
        
        function setSigmaFromOutlierDistribution(self,rejectionPDFRatio)
            % The function sets sigma (which is the initial seed in the
            % IRLS algorithm) based on the empirically determined pdf for
            % the outliers. A rejectionPDFRatio = 1e5 means that the
            % outlier pdf to noise pdf ratio is 1e5, and only point further
            % out than that are given a different sigma.
            %
            % This works remarkably well for the real GPS data, but only a
            % very small effect on the synthetic data. We speculate this is
            % because the outlier errors are correlated in the GPS data.
            if nargin < 2
                rejectionPDFRatio = 1e5;
            end
            self.setToFullTensionWithIteratedIQAD();

            [self.outlierDistribution,self.alpha] = RobustSmoothingSpline.estimateOutlierDistributionFromKnownNoise(reshape(self.epsilon,[],1),self.noiseDistribution);
            self.noiseDistanceDistribution = TwoDimDistanceDistribution(self.noiseDistribution);
            self.outlierDistanceDistribution = TwoDimDistanceDistribution(self.outlierDistribution);
            if self.alpha > 0
                self.sigma = RobustSmoothingSpline.sigmaFromOutlierDistribution(self.alpha,self.outlierDistanceDistribution,self.noiseDistanceDistribution,self.epsilon_d,rejectionPDFRatio);
                self.spline_x.sigma = self.sigma;
                self.spline_y.sigma = self.sigma;
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Minimization
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [lambda,fval] = minimize(self,penaltyFunction)
            % the penalty function *must* take a bivariate tension spline
            % object and return a scalar. The function will be minimized by
            % varying lambda.
            [~,lambdaBelow,lambdaAbove] = SmoothingSpline.minimizeFunctionOfSplineWithGridSearch(self,penaltyFunction);
            lambda = SmoothingSpline.minimizeFunctionOfSplineBounded(self,penaltyFunction,lambdaBelow,lambdaAbove);
            fval = penaltyFunction(self);
        end
        
        function [lambda,mse] = minimizeExpectedMeanSquareErrorInPercentileRange(self,pctmax)
            zmin = 0;
            zmax = self.noiseDistanceDistribution.locationOfCDFPercentile(pctmax);
            expectedVariance = self.noiseDistanceDistribution.varianceInRange(zmin,zmax);
            [lambda,mse] = self.minimize( @(aSmoothingSpline) aSmoothingSpline.expectedMeanSquareErrorInDistanceRange(zmin,zmax,expectedVariance) );
        end
        
        function [lambda,mse] = minimizeExpectedMeanSquareErrorInNoiseRange(self)
            if ~isempty(self.alpha) && self.alpha > 0
                [zoutlier,expectedVariance] = RobustSmoothingSpline.locationOfNoiseToOutlierPDFRatio(self.alpha,self.outlierDistanceDistribution,self.noiseDistanceDistribution,1);
                [lambda,mse] =self.minimize( @(spline) spline.expectedMeanSquareErrorInDistanceRange(0,zoutlier,expectedVariance));
            else
                [lambda,mse] =self.minimize( @(spline) spline.expectedMeanSquareError );
            end
        end
        
        
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Superclass override
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        function varargout = subsref(self, index)
            %% Subscript overload
            %
            % The forces subscript notation to behave as if it is
            % evaluating a function.
            idx = index(1).subs;
            switch index(1).type
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FEVAL / COMPOSE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                case '()'
                    if length(idx) >= 1
                        time = idx{1};
                    end
                    
                    if length(idx) >= 2
                        NumDerivatives = idx{2};
                    else
                        NumDerivatives = 0;
                    end
                    
                    varargout{1} = self.spline_x(time,NumDerivatives) + self.spline_mean_x(time,NumDerivatives);
                    varargout{2} = self.spline_y(time,NumDerivatives) + self.spline_mean_y(time,NumDerivatives);
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GET %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                case '.'
                    [varargout{1:nargout}] = builtin('subsref',self,index);
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RESTRICT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                case '{}'
                    error('The BivariateSmoothingSpline class does not know what to do with {}.');
                otherwise
                    error('Unexpected syntax');
            end
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Methods for setting *tension* (x_T) instead of lambda.
        %
        % If the process is isotropic, tension and lambda should be same.
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [lambda_x,lambda_y] = lambdasForTension(self,x_T)
            lambda_x = NaN; lambda_y = NaN;
            if isempty(self.tensionLambdaTable)
                return;
            else
                index = find( self.tensionLambdaTable(:,1) > x_T/1.1 & self.tensionLambdaTable(:,1) < x_T*1.1,1,'first' );
                if isempty(index)
                    return;
                else
                    lambda_x = self.tensionLambdaTable(index,2);
                    lambda_y = self.tensionLambdaTable(index,3);
                end
            end
        end
        
        function insertLambdasForTension(self,x_T,lambda_x,lambda_y)
            self.tensionLambdaTable = cat(1,self.tensionLambdaTable,[x_T, lambda_x, lambda_y]);
            [~,i] = sort(self.tensionLambdaTable(:,1));
            self.tensionLambdaTable = self.tensionLambdaTable(i,:);
        end
        
        function set.tensionValue(self,x_T)
            [lambda_x,lambda_y] = self.lambdasForTension(x_T);
            if any(isnan([lambda_x,lambda_y]))
                self.spline_x.tensionValue = x_T;
                self.spline_y.tensionValue = x_T;
                self.insertLambdasForTension(x_T, self.spline_x.lambda, self.spline_y.lambda);
            else
                self.spline_x.lambda = lambda_x;
                self.spline_y.lambda = lambda_y;
            end
        end
        
        function x_T = get.tensionValue(self)
            x_T = std(self.spline_x.uniqueValuesAtHighestDerivative());
        end
        
        function tensionValue = minimizeTension(self,penaltyFunction)
            % the penalty function *must* take a bivariate tension spline
            % object and return a scalar. The function will be minimized by
            % varying x_T across both directions
            [~,tensionValueBelow,tensionValueAbove] = GPSSmoothingSpline.minimizeFunctionOfSplineWithGridSearch(self,penaltyFunction);
            tensionValue = mean([tensionValueAbove,tensionValueBelow]);
            self.tensionValue = tensionValue;
            %             tensionValue = GPSSmoothingSpline.minimizeFunctionOfSplineBounded(self,penaltyFunction,tensionValueBelow,tensionValueAbove);
        end
        
    end
    
    methods (Static)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Minimization on the *tension* directly (not lambda)
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [lambda,lambdaBelow,lambdaAbove] = minimizeFunctionOfBivariateSplineWithGridSearch(aBivariateSpline,functionOfBivariateSpline,dlog10tensionValue,n)
            tensionValues = aBivariateSpline.tensionValue;
            err = functionOfBivariateSpline(aBivariateSpline);
            isConstrained = aBivariateSpline.isConstrained;
            nLambdas = 1;
            index = 1;
            
            if nargin < 4
                dlog10tensionValue = 0.5; % step size in factors of 10
                n = 3; % number of steps to make each direction
            end
            
            while ( (index <= 2 && isConstrained(index) ~=1) || index >= nLambdas-1 )   
                if index <= 2 && isConstrained(index) ~=1
                    % expand search below
                    tensionValue0 = tensionValues(1);
                    tensionValues = cat(1,10.^linspace(log10(tensionValue0)-n*dlog10tensionValue,log10(tensionValue0)-dlog10tensionValue,n)',tensionValues);
                    err = cat(1,zeros(n,1),err);
                    isConstrained = cat(1,zeros(n,1),isConstrained);
                    for iLambda = 1:n
                        aBivariateSpline.tensionValue = tensionValues(iLambda);
                        err(iLambda) =  functionOfBivariateSpline(aBivariateSpline);
                        isConstrained(iLambda) = aBivariateSpline.isConstrained;
                    end
                    
                    nLambdas = length(tensionValues);
                    index = index+n;
                end
                
                if index >= nLambdas-1
                    % expand search above
                    tensionValue0 = tensionValues(end);
                    tensionValues = cat(1,tensionValues,10.^linspace(log10(tensionValue0)+dlog10tensionValue,log10(tensionValue0)+n*dlog10tensionValue,n)');
                    err = cat(1,err,zeros(n,1));
                    isConstrained = cat(1,isConstrained,zeros(n,1));
                    for iLambda = (nLambdas+1):(nLambdas+n)
                        aBivariateSpline.tensionValue = tensionValues(iLambda);
                        err(iLambda) =  functionOfBivariateSpline(aBivariateSpline);
                        isConstrained(iLambda) = aBivariateSpline.isConstrained;
                    end
                    
                    nLambdas = length(tensionValues);
                end
                
                [~,index] = min(err);
            end
            lambda = tensionValues(index);
            aBivariateSpline.tensionValue = lambda; 
            lambdaBelow = tensionValues(index-1);
            if length(lambda)>index
                lambdaAbove = tensionValues(index+1);
            else
                lambdaAbove = tensionValues(index);
            end
        end
        
        function tensionValue = minimizeFunctionOfSplineBounded(aBivariateSpline,functionOfBivariateSpline,x1,x2)
            epsilon = 1e-15;
            errorFunction = @(log10tensionValuePlusEpsilon) GPSSmoothingSpline.functionOfSplineWrapper(aBivariateSpline,log10tensionValuePlusEpsilon,functionOfBivariateSpline);
            optimalLog10tensionValuePlusEpsilon = fminbnd( errorFunction, log10(x1+epsilon), log10(x2+epsilon), optimset('TolX', 0.01, 'TolFun', 0.001) );
            tensionValue = 10^optimalLog10tensionValuePlusEpsilon - epsilon;
            aBivariateSpline.tensionValue = tensionValue;
        end
        
        function error = functionOfSplineWrapper(aBivariateSpline, log10tensionValuePlusEpsilon, functionOfSpline)
            epsilon = 1e-15;
            aBivariateSpline.tensionValue = 10^log10tensionValuePlusEpsilon-epsilon;
            error = functionOfSpline(aBivariateSpline);
        end
        
    end
    
end