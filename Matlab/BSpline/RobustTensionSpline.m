classdef RobustTensionSpline < TensionSpline
    %RobustTensionSpline Allows for outliers
    
    properties
        % These are not used directly in minimization, only the
        % 'distribution' property in the TensionSpline superclass is used
        % directly.
        noiseDistribution
        outlierDistribution
        alpha
        
        w_epsilon   % weighting of the errors
        
        zmin
        zmax
    end
    
    properties (Dependent)
       t_all
       x_all
    end
    
    methods
        function self = RobustTensionSpline(t,x,distribution,varargin)
            % Construct a new 'robust' distribution by adding a Student's
            % t-distribution with 1000 times the variance, but assuming
            % only .01 percent outliers.
            % Override any user settings for lambda
            didOverrideLambda = 0;
            alpha = 1/10000;
            outlierMethod = OutlierMethod.sigmaMethod;
            rejectionPDFRatio = 1e5;
            minimizationPDFRatio = 1;
            lambdaArgument = Lambda.optimalIterated;
            for k = 1:2:length(varargin)
                if strcmp(varargin{k}, 'lambda')
                    lambdaArgument = varargin{k+1};
                    varargin{k+1} = Lambda.fullTensionExpected;
                    didOverrideLambda = 1;
                elseif strcmp(varargin{k}, 'alpha')
                    alpha = varargin{k+1};
                elseif strcmp(varargin{k}, 'outlierMethod')
                    outlierMethod = varargin{k+1};
                elseif strcmp(varargin{k}, 'rejectionPDFRatio')
                    rejectionPDFRatio = varargin{k+1};
                elseif strcmp(varargin{k}, 'minimizationPDFRatio')
                    minimizationPDFRatio = varargin{k+1};
                end
            end
            if didOverrideLambda == 0
                varargin{end+1} = 'lambda';
                varargin{end+1} = Lambda.fullTensionExpected;
            end
            
            noiseDistribution = distribution;
            nu = 3.0;
            sigma = sqrt(noiseDistribution.variance*1000*(nu-2)/nu);
            outlierDistribution = StudentTDistribution(sigma,nu);
            distribution = AddedDistribution(alpha,outlierDistribution,noiseDistribution);
            
            self@TensionSpline(t,x,distribution,varargin{:});
            
            self.noiseDistribution = noiseDistribution;
            self.outlierDistribution = outlierDistribution;
            self.alpha = alpha;
            self.outlierThreshold = self.noiseDistribution.locationOfCDFPercentile(1-1/10000/2);
            self.w_epsilon = ones(size(self.t));
            
            switch outlierMethod
                case OutlierMethod.sigmaMethod
                    self.setSigmaFromOutlierDistribution(rejectionPDFRatio);
                case OutlierMethod.distributionMethod
                    self.setDistributionFromOutlierDistribution();
                case OutlierMethod.distributionWeightAndSigmaMethod
                    self.setDistributionWeightAndSigmaFromOutlierDistribution(rejectionPDFRatio);
                case OutlierMethod.knotRemovalMethod
                    self.removeOutlierKnotsUsingOutlierDistribution(rejectionPDFRatio);
                case OutlierMethod.weightingMethod
                    self.generateEpsilonWeightingFromOutlierDistribution();
            end
                            
            if isenum(lambdaArgument)
                switch lambdaArgument
                    case {Lambda.optimalExpected,Lambda.fullTensionExpected}
                        fprintf('Ignoring your request for expected tension values. Tension set to fullTensionIterated.');
                    case Lambda.fullTensionIterated
                        
                    case {Lambda.optimalIterated}
                        if outlierMethod == OutlierMethod.weightingMethod
                            self.minimize( @(spline) spline.expectedMeanSquareErrorWithWeighting(self.noiseDistribution.variance) );
                        else
                            self.minimizeExpectedMeanSquareErrorInNoiseRange(minimizationPDFRatio);
                        end
                end
            elseif isscalar(lambdaArgument)
                self.lambda = lambdaArgument;
            else
                error('Invalid choice for lambda. Lambda must be either a scalar or the enum Lambda.');
            end
                
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Overriding the expected MSE minimization algorithm
        %
        % 2019/02/26 - 
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         function minimizeExpectedMeanSquareError(self)
%             
%             % This value of beta is chosen from
%             % MakeTableOutliersWithAddedDistributionBlind.m as a good value
%             % for distribution with very few outliers.
%             beta = 1/400;
%             zmin_ = self.noiseDistribution.locationOfCDFPercentile(beta/2);
%             zmax_ = self.noiseDistribution.locationOfCDFPercentile(1-beta/2);
%             
%             if ~(isempty(self.outlierDistribution) || self.alpha < 0.01)
%                 % If we a non-trivial outlier distribution, includepoints
%                 % that are more likely than not to be part of the
%                 % characterized noise distribution.
%                 noiseOdds = 3; % This is from MakeTableOutliersBlindMinimization.m
%                 f = @(z) abs( (1-self.alpha)*self.noiseDistribution.pdf(z) - noiseOdds*self.alpha*self.outlierDistribution.pdf(z) );
%                 zoutlier = abs(fminsearch(f,sqrt(self.noiseDistribution.variance)));
%                 if zoutlier < zmax_
%                     zmin_ = -zoutlier;
%                     zmax_ = zoutlier;
%                 end
%             end
%             self.minimize( @(spline) spline.expectedMeanSquareErrorInRange(zmin_,zmax_) );
%             %             fprintf('Minimizing in range: (%.1f, %.1f) with expected variance %.1f m^2\n',zmin_,zmax_,self.distribution.varianceInRange(zmin_,zmax_));
%         end
        
        function minimizeExpectedMeanSquareErrorInNoiseRange(self,minimizationPDFRatio)
           if nargin < 2
               minimizationPDFRatio = 1;
           end
           if self.alpha > 0
               [zoutlier,expectedVariance] = RobustTensionSpline.locationOfNoiseToOutlierPDFRatio(self.alpha,self.outlierDistribution,self.noiseDistribution,minimizationPDFRatio);
               self.minimize( @(spline) spline.expectedMeanSquareErrorInRange(-zoutlier,zoutlier,expectedVariance));
           else
               self.minimizeExpectedMeanSquareError();
           end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Methods for achieve full tension
        %
        % 2019/02/26 - The preferred method for now is the iterated
        % interquartile Anderson-Darling test.
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function setToFullTensionWithIteratedIQAD(self)
            % Set to full tension by minimizing the Anderson-Darling (AD)
            % error on the interquartile (IQ) range of epsilons. At each
            % iteration the outlier distribution is estimated and used to
            % refine the tension.
            self.minimize(@(spline) self.noiseDistribution.andersonDarlingInterquartileError(spline.epsilon));
            lastAlpha = 0.0;
            totalIterations = 0;
            [newOutlierDistribution, newAlpha] = RobustTensionSpline.estimateOutlierDistributionFromKnownNoise(self.epsilon,self.noiseDistribution);
            while (abs(lastAlpha-newAlpha) > 0.01 && totalIterations < 10)
                if newAlpha > 0 && ~isempty(newOutlierDistribution)
                    addedDistribution = AddedDistribution(newAlpha,newOutlierDistribution,self.noiseDistribution);
                else
                    addedDistribution = self.noiseDistribution;
                end
                self.minimize(@(spline) addedDistribution.andersonDarlingInterquartileError(spline.epsilon));
                lastAlpha = newAlpha;
                [newOutlierDistribution, newAlpha] = RobustTensionSpline.estimateOutlierDistributionFromKnownNoise(self.epsilon,self.noiseDistribution);
                totalIterations = totalIterations + 1;
            end
        end
        
        function setToFullTensionWithIteratedIQSV(self,extraVarianceFactor)
            % Set to full tension by matching the expected sample variance
            % (SV) with the interquartile (IQ) sample variance. At each
            % iteration the outlier distribution is estimated and used to
            % refine the tension.
            if nargin < 2
                extraVarianceFactor = 1.0;
            end
            pctRange = 1/2; % alpha=1/2 is the inner half of the distribution
            self.setInnerVarianceToExpectedValue(pctRange,extraVarianceFactor*self.noiseDistribution.varianceInPercentileRange(pctRange/2,1-pctRange/2));
            lastAlpha = 0.0;
            totalIterations = 0;
            [newOutlierDistribution, newAlpha] = RobustTensionSpline.estimateOutlierDistributionFromKnownNoise(self.epsilon,self.noiseDistribution);
            while (abs(lastAlpha-newAlpha) > 0.01 && totalIterations < 10)
                addedDistribution = AddedDistribution(newAlpha,newOutlierDistribution,self.noiseDistribution);
                self.setInnerVarianceToExpectedValue(pctRange,extraVarianceFactor*addedDistribution.varianceInPercentileRange(pctRange/2,1-pctRange/2));
                lastAlpha = newAlpha;
                [newOutlierDistribution, newAlpha] = RobustTensionSpline.estimateOutlierDistributionFromKnownNoise(self.epsilon,self.noiseDistribution);
                totalIterations = totalIterations + 1;
            end
        end
        
        function setToFullTensionWithKS(self,alpha)
            % Set to full tension by minimizing the Kolmogorov-Smirnof
            % error on epsilons found to be within the percentile range set
            % by alpha.
            if nargin == 1
                alpha = 1/10;
            end
            
            epsilon_min = self.noiseDistribution.locationOfCDFPercentile(alpha/2);
            epsilon_max = self.noiseDistribution.locationOfCDFPercentile(1-alpha/2);
            self.minimize( @(spline) self.noiseDistribution.kolmogorovSmirnovError(spline.epsilon,epsilon_min,epsilon_max) );
        end
        
        function setToFullTensionWithSV(self,alpha)
            % This gives us "full tension" uby matching the sample variance
            % to expected variance of the noise within some percentile
            % range. Alpha sets the percentile range (e.g., 1/2 ignores any
            % points outside the distances set by alpha).
            %
            % Note that this method is *not* interquartile--the number of
            % points being used is not fixed and could even go to zero.
            if nargin == 1
                alpha = 1/10;
            end

            zmin_ = self.noiseDistribution.locationOfCDFPercentile(alpha/2);
            zmax_ = self.noiseDistribution.locationOfCDFPercentile(1-alpha/2);
            expectedVarianceInRange = self.noiseDistribution.varianceInRange(zmin_,zmax_);
            self.minimize(@(spline) abs(spline.sampleVarianceInRange(zmin_,zmax_) - expectedVarianceInRange));
        end
        
        function setToFullTensionWithInnerSVOnNoiseDistribution(self,alpha)
            % This algorithm takes the interquartile range, and requires it
            % match the expected variance in that range. This is different than
            % above in that it has a fixed number of points.
            
            zmin_ = self.noiseDistribution.locationOfCDFPercentile(alpha/2);
            zmax_ = self.noiseDistribution.locationOfCDFPercentile(1-alpha/2);
            expectedVarianceInRange = self.noiseDistribution.varianceInRange(zmin_,zmax_);
            self.minimize(@(spline) abs(spline.sampleVarianceIQ(alpha) - expectedVarianceInRange));
        end
        
        function setToFullTensionWithIQSV(self,alpha)
            % This algorithm takes the interquartile range, and requires it
            % match the expected variance in that range. This is different than
            % above in that it has a fixed number of points.
            
            zmin_ = self.distribution.locationOfCDFPercentile(alpha/2);
            zmax_ = self.distribution.locationOfCDFPercentile(1-alpha/2);
            expectedVarianceInRange = self.distribution.varianceInRange(zmin_,zmax_);
            self.setInterquartileSampleVarianceToExpectedValue(alpha,expectedVarianceInRange);
        end
        
        function setInterquartileSampleVarianceToExpectedValue(self,alpha,expectedVarianceInRange)
            self.minimize(@(spline) abs(spline.sampleVarianceIQ(alpha) - expectedVarianceInRange));
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
            epsilon_full = self.epsilon;
            [self.outlierDistribution,self.alpha] = RobustTensionSpline.estimateOutlierDistributionFromKnownNoise(epsilon_full,self.noiseDistribution);
            if self.alpha > 0
                self.sigma = RobustTensionSpline.sigmaFromOutlierDistribution(self.alpha,self.outlierDistribution,self.noiseDistribution,epsilon_full,rejectionPDFRatio);
            end
        end
        
        function setDistributionFromOutlierDistribution(self)
            % This function set the distribution property (which is used
            % for minimization) from the outlier distribution.
            %
            % This strategy works remarkably badly, even with the correct
            % outlier distribution.
            self.setToFullTensionWithIteratedIQAD();
            epsilon_full = self.epsilon;
            [self.outlierDistribution,self.alpha] = RobustTensionSpline.estimateOutlierDistributionFromKnownNoise(epsilon_full,self.noiseDistribution);  
            if self.alpha > 0.0
                self.distribution = AddedDistribution(self.alpha,self.outlierDistribution,self.noiseDistribution);
            end
        end
        
        function setDistributionWeightAndSigmaFromOutlierDistribution(self,rejectionPDFRatio)
            % This function set the distribution property (which is used
            % for minimization) from the outlier distribution. It also sets
            % the weight function and sigma to match the outlier
            % distribution for points above the rejectionPDFRatio.
            %
            % This strategy works okay... but far more variable than just
            % setting sigma
            if nargin < 2
                rejectionPDFRatio = 1e5;
            end
            self.setToFullTensionWithIteratedIQAD();
            epsilon_full = self.epsilon;
            [self.outlierDistribution,self.alpha] = RobustTensionSpline.estimateOutlierDistributionFromKnownNoise(epsilon_full,self.noiseDistribution);            
            if self.alpha > 0.0
                self.distribution = AddedDistribution(self.alpha,self.outlierDistribution,self.noiseDistribution);
                f = @(z) abs( (self.alpha/(1-self.alpha))*self.outlierDistribution.pdf(-abs(z))/self.noiseDistribution.pdf(-abs(z)) - rejectionPDFRatio);
                z_outlier = fminsearch(f,sqrt(self.noiseDistribution.variance));
                noiseIndices = abs(epsilon_full) <= z_outlier;
                
                self.distribution.w = @(z) noiseIndices .* self.noiseDistribution.w(z) + (~noiseIndices) .* self.outlierDistribution.w(z);
                self.sigma = noiseIndices .* sqrt(self.noiseDistribution.variance) + (~noiseIndices) .* sqrt(self.outlierDistribution.variance);
            end
        end
        
        function removeOutlierKnotsUsingOutlierDistribution(self,rejectionPDFRatio)
            % The strategy here is to remove support from points deemed to
            % be outliers. So we set to full tension, identify outliers
            % (using the default threshold), and then retension
            self.setToFullTensionWithIteratedIQAD();
            epsilon_full = self.epsilon;
            [self.outlierDistribution,self.alpha] = RobustTensionSpline.estimateOutlierDistributionFromKnownNoise(epsilon_full,self.noiseDistribution);
            
            if self.alpha < 0
                return;
            end

            f = @(z) abs( (self.alpha/(1-self.alpha))*self.outlierDistribution.pdf(-abs(z))/self.noiseDistribution.pdf(-abs(z)) - rejectionPDFRatio);
            z_outlier = fminsearch(f,sqrt(self.noiseDistribution.variance));
            newKnotIndices = find(abs(epsilon_full) <= z_outlier);

            if newKnotIndices(1) ~= 1
                newKnotIndices = cat(2,1,newKnotIndices);
            end
            if newKnotIndices(end) ~= length(self.t)
                newKnotIndices = cat(2,newKnotIndices,length(self.t));
            end
            
            self.t_knot = InterpolatingSpline.KnotPointsForPoints(self.t(newKnotIndices),self.K,self.knot_dof);
            self.X = BSpline.Spline( self.t, self.t_knot, self.K, 0 ); % NxM
            
            % Now we need a quadrature (integration) grid that is finer
            % if S=T we can optimize this much better because it's all
            % piecewise constant.
            Q = 10*length(self.t); % number of points on the quadrature grid
            tq = linspace(self.t(1),self.t(end),Q)';
            B = BSpline.Spline( tq, self.t_knot, self.K, self.T );
            self.V = squeeze(B(:,:,self.T+1)); % QxM
            
            % Precompute some matrices that might be used again later,
            [self.XWX,self.XWx,self.VV] = TensionSpline.PrecomputeTensionSolutionMatrices(self.X,self.V,sqrt(self.distribution.variance),self.x);
            
            self.tensionParameterDidChange();
        end
        
        function generateEpsilonWeightingFromOutlierDistribution(self)
            % This is a fairly experiment method that works okay, not
            % great, but good enough to keep for the record. The idea is
            % computed the expected mean square error *weighted* with the
            % probability that a given point is part of the expected
            % distribution.
            %
            % This method builds the weighting function, and the
            % expectedMeanSquareErrorWithWeighting function does the actual
            % minimization.
            self.setToFullTensionWithIteratedIQAD();
            epsilon_full = self.epsilon;
            [self.outlierDistribution,self.alpha] = RobustTensionSpline.estimateOutlierDistributionFromKnownNoise(epsilon_full,self.noiseDistribution);
            if self.alpha > 0.0
                self.distribution = AddedDistribution(self.alpha,self.outlierDistribution,self.noiseDistribution);
                self.w_epsilon = (1-self.alpha)*self.noiseDistribution.pdf(epsilon_full)./self.distribution.pdf(epsilon_full);
            else
                self.w_epsilon = ones(size(self.t));
            end
        end
        
        function [MSE, n] = expectedMeanSquareErrorWithWeighting(self,expectedVariance)
            epsilon = self.epsilon;
            n = sum(self.w_epsilon);
            X2 = sum((epsilon).^2 .* self.w_epsilon)/n;
            
            Sii = diag(self.smoothingMatrix);
            traceS = sum( Sii .* self.w_epsilon );
            
            MSE = X2/expectedVariance + 2*traceS/n - 1;
        end
        
    end
    
    methods (Static)
        
        function sigma = sigmaFromOutlierDistribution(alpha,outlierDistribution,noiseDistribution,epsilon_full,outlierOdds)
            % Returns an array that serves as the initial seed in the IRLS
            % algorithm. This *usually* doesn't matter, but does have an
            % effect in some case (when the outliers are correlated?).
            f = @(z) abs( (alpha/(1-alpha))*outlierDistribution.pdf(-abs(z))/noiseDistribution.pdf(-abs(z)) - outlierOdds);
            z_outlier = fminsearch(f,sqrt(noiseDistribution.variance));
            noiseIndices = abs(epsilon_full) <= z_outlier;           
            sigma = noiseIndices .* sqrt(noiseDistribution.variance) + (~noiseIndices) .* sqrt(outlierDistribution.variance);
        end
        
        function [z,expectedVariance] = varianceCrossOverFromOutlierDistribution(alpha,outlierDistribution,noiseDistribution)
            % Given the same *total* variance, where (in z) do two
            % distributions with different ratios of outliers have the same
            % integrated variance?
            %
            % This appears to serve as a reasonable cutoff for minimizing
            % the expected mean square error.
            addedDistribution = AddedDistribution(alpha,outlierDistribution,noiseDistribution);
            weakAddedDistribution = AddedDistribution(alpha/1.5,StudentTDistribution(sqrt((addedDistribution.variance-(1-alpha/1.5)*noiseDistribution.variance)/(3*alpha/2)),3.0),noiseDistribution);
            strongAddedDistribution = AddedDistribution(1.5*alpha,StudentTDistribution(sqrt((addedDistribution.variance-(1-alpha*1.5)*noiseDistribution.variance)/(3*alpha*1.5)),3.0),noiseDistribution);
            f = @(z) abs(weakAddedDistribution.varianceInRange(0,z)-strongAddedDistribution.varianceInRange(0,z));
            z = fminbnd(f,sqrt(noiseDistribution.variance),sqrt(strongAddedDistribution.variance));
            expectedVariance = addedDistribution.varianceInRange(-z,z);
        end
        
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
        
    end
end

