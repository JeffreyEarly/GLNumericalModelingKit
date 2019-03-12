classdef BivariateTensionSpline < handle
    % BivariateTensionSpline Fit 2D (two-dependent, one-independent) data
    % isotropically.
    %
    %   4 argument initialization
    %       f = BivariateTensionSpline(t,x,y,distribution);
    %   where
    %       t               time (seconds or datetime)
    %       x               latitude
    %       y               longitude
    %       distribution    distribution of errors
    %       f               spline interpolant
    %
    %   BivariateTensionSpline takes a number of optional input argument pairs.
    
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
        % 'distribution' property in the TensionSpline superclass is used
        % directly.
        noiseDistribution
        outlierDistribution
        alpha
        noiseDistanceDistribution
        outlierDistanceDistribution
        sigma
        
        % splines
        spline_constrained_x
        spline_constrained_y
        spline_x
        spline_y
    end
    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Initialization
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function self = BivariateTensionSpline(t,x,y,distribution,varargin)
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
            
            t_knot = cat(1,min(t)*ones(self.K+1,1),max(t)*ones(self.K+1,1));
            self.spline_constrained_x = ConstrainedSpline(self.t,self.x,self.K+1,t_knot,self.distribution,[]);
            self.spline_constrained_y = ConstrainedSpline(self.t,self.y,self.K+1,t_knot,self.distribution,[]);
            
            self.x_bar = self.spline_constrained_x(self.t);
            self.y_bar = self.spline_constrained_y(self.t);
            
            self.x_prime = self.x - self.x_bar;
            self.y_prime = self.y - self.y_bar;
             
            if self.shouldUseRobustFit == 1
                self.spline_x = RobustTensionSpline(self.t,self.x_prime,self.distribution,'K',self.K,'T',self.T,'lambda',Lambda.fullTensionIterated);
                self.spline_y = RobustTensionSpline(self.t,self.y_prime,self.distribution,'K',self.K,'T',self.T,'lambda',self.spline_x.lambda);
            else
                self.spline_x = TensionSpline(self.t,self.x_prime,self.distribution,'K',self.K,'T',self.T);
                self.spline_y = TensionSpline(self.t,self.y_prime,self.distribution,'K',self.K,'T',self.T,'lambda',self.spline_x.lambda);
            end
            self.lambda = self.spline_x.lambda;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Evaluation
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [x,y] = xyAtTime(self,time)
            x = self.spline_x(time) + self.spline_constrained_x(time);
            y = self.spline_y(time) + self.spline_constrained_y(time);
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
        % Measures of expected error
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function MSE = expectedMeanSquareError(self)
            % This is the *expected* mean-square error normalized by the
            % variance, simply added together across dimensions. Does *not*
            % include the mean.
            
            MSE = self.spline_x.MSE + self.spline_y.MSE;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Measures of expected error?restricted to a subset of points
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
            indices = find(epsilon >= zmin & epsilon <= zmax);
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
            self.minimize(@(spline) self.distribution.andersonDarlingInterquartileError( reshape(spline.epsilon,[],1) ));
            
            lastAlpha = 0.0;
            totalIterations = 0;
            [newOutlierDistribution, newAlpha] = RobustTensionSpline.estimateOutlierDistributionFromKnownNoise(reshape(self.epsilon,[],1),self.distribution);
            while (abs(lastAlpha-newAlpha) > 0.01 && totalIterations < 10)
                if newAlpha > 0 && ~isempty(newOutlierDistribution)
                    addedDistribution = AddedDistribution(newAlpha,newOutlierDistribution,self.distribution);
                else
                    addedDistribution = self.distribution;
                end
                self.minimize(@(spline) addedDistribution.andersonDarlingInterquartileError( reshape(spline.epsilon,[],1) ));
                lastAlpha = newAlpha;
                [newOutlierDistribution, newAlpha] = RobustTensionSpline.estimateOutlierDistributionFromKnownNoise(reshape(self.epsilon,[],1),self.distribution);
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

            [self.outlierDistribution,self.alpha] = RobustTensionSpline.estimateOutlierDistributionFromKnownNoise(reshape(self.epsilon,[],1),self.noiseDistribution);
            self.noiseDistanceDistribution = TwoDimDistanceDistribution(self.noiseDistribution);
            self.outlierDistanceDistribution = TwoDimDistanceDistribution(self.outlierDistribution);
            if self.alpha > 0
                self.sigma = RobustTensionSpline.sigmaFromOutlierDistribution(self.alpha,self.outlierDistanceDistribution,self.noiseDistanceDistribution,self.epsilon_d,rejectionPDFRatio);
                self.spline_x.sigma = self.sigma;
                self.spline_y.sigma = self.sigma;
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Minimization
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function lambda = minimize(self,penaltyFunction)
            % the penalty function *must* take a bivariate tension spline
            % object and return a scalar. The function will be minimized by
            % varying lambda.
            [~,lambdaBelow,lambdaAbove] = TensionSpline.minimizeFunctionOfSplineWithGridSearch(self,penaltyFunction);
            lambda = TensionSpline.minimizeFunctionOfSplineBounded(self,penaltyFunction,lambdaBelow,lambdaAbove);
        end
        
        function minimizeExpectedMeanSquareErrorInNoiseRange(self,minimizationPDFRatio)
            if nargin < 2
                minimizationPDFRatio = 1;
            end
            if self.alpha > 0
                [zoutlier,expectedVariance] = RobustTensionSpline.locationOfNoiseToOutlierPDFRatio(self.alpha,self.outlierDistanceDistribution,self.noiseDistanceDistribution,minimizationPDFRatio);
                expectedVariance2D = AddedDistribution(self.alpha,self.noiseDistribution,self.outlierDistribution).varianceInRange(-zoutlier,zoutlier);
                fprintf('Expected variance (1d/2d): %.1f / %.1f \n',expectedVariance,expectedVariance2D);
                self.minimize( @(spline) spline.expectedMeanSquareErrorInDistanceRange(0,zoutlier,expectedVariance));
            else
                self.minimizeExpectedMeanSquareError();
            end
        end
        
    end
    
end