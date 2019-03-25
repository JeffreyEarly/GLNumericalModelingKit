classdef RobustTensionSpline < TensionSpline
    %RobustTensionSpline Allows for outliers to occur.
    % 
    
    properties
        % These are not used directly in minimization, only the
        % 'distribution' property in the TensionSpline superclass is used
        % directly.
        noiseDistribution
        
        w_epsilon   % weighting of the errors
    end

    methods
        function self = RobustTensionSpline(t,x,distribution,varargin)
            % Construct a new 'robust' distribution by adding a Student's
            % t-distribution with 1000 times the variance, but assuming
            % only .01 percent outliers.
            % Override any user settings for lambda
            didOverrideLambda = 0;
            alpha = 1/10000;
            outlierMethod = OutlierMethod.sigmaFullTensionMethod;
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
            if alpha > 0
                distribution = AddedDistribution(alpha,outlierDistribution,noiseDistribution);
            end
            
            self@TensionSpline(t,x,distribution,varargin{:});
            
            self.noiseDistribution = noiseDistribution;
            self.outlierDistribution = outlierDistribution;
            self.alpha = alpha;
            self.outlierThreshold = self.noiseDistribution.locationOfCDFPercentile(1-1/10000/2);
            self.w_epsilon = ones(size(self.t));
            
            
            % Default algorithm is as follows:
            % 1. estimate the outlier distribution by going to full tension
            % 2. find the ranged expected mean square error (using the pdf
            %    cross over point)
            % 3. apply some user selected outlier method
            % 4. re-compute the ranged expected mean square error.
            % 5. choose the global minimum (from step 2 and 4)
            
            
            % Step 1---estimate the outlier distribution *always*
            [self.outlierDistribution, self.alpha] = self.setToFullTensionWithIteratedIQAD();
            self.lambdaAtFullTension = self.lambda;
            self.sigmaAtFullTension = self.distribution.w(self.epsilon);
            
            % Step 2---minimize under current conditions
        
            if isenum(lambdaArgument)
                if lambdaArgument == Lambda.optimalIterated || lambdaArgument == Lambda.optimalExpected
                    if outlierMethod ~= OutlierMethod.weightingMethod
                        mse1 = self.minimizeExpectedMeanSquareErrorInNoiseRange(minimizationPDFRatio);
                        mse1_lambda = self.lambda;
                    end
                end
            end
            
            % Step 3---apply outlier method
            switch outlierMethod
                case OutlierMethod.sigmaFullTensionMethod
                    self.sigma = self.sigmaAtFullTension;
                case OutlierMethod.sigmaRejectionMethod
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
            
            % Step 4---re-compute minimum expected mean square error
            if isenum(lambdaArgument)
                switch lambdaArgument
                    case {Lambda.fullTensionIterated,Lambda.fullTensionExpected}
                        self.lambda = self.lambdaAtFullTension;
                    case {Lambda.optimalIterated,Lambda.optimalExpected}
                        if outlierMethod == OutlierMethod.weightingMethod
                            self.minimize( @(spline) spline.expectedMeanSquareErrorWithWeighting(self.noiseDistribution.variance) );
                        else
                            mse2 = self.minimizeExpectedMeanSquareErrorInNoiseRange(minimizationPDFRatio);
                            
                            % Step 5---set to global minimum
                            if mse1 < mse2
                                self.sigma = sqrt(self.distribution.variance);
                                self.lambda = mse1_lambda;
                            end
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
        
        function mse = minimizeExpectedMeanSquareErrorInNoiseRange(self,minimizationPDFRatio)
           if nargin < 2
               minimizationPDFRatio = 1;
           end
           if self.alpha > 0
               [zoutlier,expectedVariance] = RobustTensionSpline.locationOfNoiseToOutlierPDFRatio(self.alpha,self.outlierDistribution,self.noiseDistribution,minimizationPDFRatio);
               self.minimize( @(spline) spline.expectedMeanSquareErrorInRange(-zoutlier,zoutlier,expectedVariance));
               mse = self.expectedMeanSquareErrorInRange(-zoutlier,zoutlier,expectedVariance);
           else
               self.minimizeExpectedMeanSquareError();
               mse = self.expectedMeanSquareError();
           end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Experimental methods
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
              
        function [MSEoutlier,MSEnoise] = findParameterizedNoiseCrossoverPoint(self)
            zmin = self.distribution.locationOfCDFPercentile(1/10000);
            zmax = self.distribution.locationOfCDFPercentile(1/10);
            z = linspace(abs(zmax),abs(zmin),10)';
            
            MSEoutlier = zeros(length(z),1);
            MSEnoise = zeros(length(z),1);
            
            epsilon = self.epsilon;
            
            for iZ = 1:length(z)
               zValue = z(iZ);
               noise_indices = find(epsilon >= -zValue & epsilon <= zValue);
               outlier_indices = setdiff(1:length(self.t),noise_indices);
               MSEoutlier(iZ) = self.expectedMeanSquareErrorFromCVForPointsAtIndices(outlier_indices);
               MSEnoise(iZ) = self.expectedMeanSquareErrorForPointsAtIndices(noise_indices,self.distribution.variance)*self.distribution.variance;
            end
        end
        
        function zCrossover = findSampleVarianceCrossover(self)            
            zCrossover = fminsearch(@(z) self.diffExpectedActualSampleVariance(z),sqrt(self.distribution.variance));
        end
        
        function dSampleVariance = diffExpectedActualSampleVariance(self,z)
            epsilon = self.epsilon;
            noise_indices = find(epsilon >= -abs(z) & epsilon <= abs(z));
            n_eff = self.effectiveSampleSizeFromVarianceOfTheMeanForIndices(noise_indices);
            dSampleVariance = abs( (1-1/n_eff)*self.distribution.variance - mean(epsilon(noise_indices).^2) );
        end
        
        function MSE = expectedMeanSquareErrorRobust(self)
            zCrossover = fminsearch(@(z) self.diffExpectedActualSampleVariance(z),sqrt(self.distribution.variance));
%             fprintf('%f\t',zCrossover);
            MSE = self.expectedMeanSquareErrorInRange(-abs(zCrossover),abs(zCrossover));
        end
        
        function [z,n_eff,expectedSampleVariance,actualSampleVariance] = findSampleVarianceCrossoverWithGridSearch(self)
            zmin = self.distribution.locationOfCDFPercentile(1/10000);
            zmax = self.distribution.locationOfCDFPercentile(1/10);
            z = linspace(abs(zmax),abs(zmin),10)';
            
            n_eff = zeros(length(z),1);
            expectedSampleVariance = zeros(length(z),1);
            actualSampleVariance = zeros(length(z),1);
            
            epsilon = self.epsilon;
            
            for iZ = 1:length(z)
                zValue = z(iZ);
                noise_indices = find(epsilon >= -zValue & epsilon <= zValue);
                
                n_eff(iZ) = self.effectiveSampleSizeFromVarianceOfTheMeanForIndices(noise_indices);
                expectedSampleVariance(iZ) = (1-1/n_eff(iZ))*self.distribution.variance;
                actualSampleVariance(iZ) = mean(epsilon(noise_indices).^2);
            end
        end
        
        function MSE = expectedMeanSquareErrorRobustOld(self,zmin,zmax)
            epsilon = self.epsilon;
            noise_indices = find(epsilon >= zmin & epsilon <= zmax);
            expectedVariance = self.distribution.varianceInRange(zmin,zmax);
            MSE = self.expectedMeanSquareErrorForPointsAtIndices(noise_indices,expectedVariance);
            
            expectedNumberOfOutliers = length(self.t)*(self.distribution.cdf(zmin) + (1-self.distribution.cdf(zmax)));
            
            outlier_indices = find(epsilon < zmin | epsilon > zmax);
            n_outliers = length(outlier_indices);
            if n_outliers > 3*expectedNumberOfOutliers
                MSE_outlier = self.expectedMeanSquareErrorFromCVForPointsAtIndices(outlier_indices);
                
                n_eff = self.effectiveSampleSizeFromVarianceOfTheMeanForIndices(outlier_indices);
                outlier_variance = mean(epsilon(outlier_indices).^2)/(1-1/n_eff);
                MSE_outlier = MSE_outlier/outlier_variance;
                
%                 alpha = n_outliers/length(self.t);
%                 MSE = (1-alpha)*MSE + alpha*MSE_outlier/outlier_variance;
                MSE = (length(noise_indices)*MSE + n_outliers*MSE_outlier)/length(self.t);
            end
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
        % Measures of expected error restricted to a subset of points
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function MSE = expectedMeanSquareErrorFromCVForPointsAtIndices(self,indices)
            % Cross-validation (CV) estimate for the mean square error from
            % Green and Silverman, equation 3.5
            epsilon = self.epsilon;
            
            S = self.smoothingMatrix;
            S = S(indices,indices);
            Sii = diag(S);
            
            MSE = mean( (epsilon(indices)./(1-Sii)).^2 );
        end
        
        function MSE = expectedMeanSquareErrorInterquartile(self)
            epsilon = self.x - self.ValueAtPoints(self.t);
            sortedEpsilon = sort(epsilon);
            n = length(epsilon);
            zmin = sortedEpsilon(floor(n/4));
            zmax = sortedEpsilon(ceil(3*n/4));
            expectedVariance = self.distribution.varianceInPercentileRange(0.25,0.75);
            MSE = self.expectedMeanSquareErrorForPointsAtIndices(epsilon >= zmin & epsilon <= zmax,expectedVariance);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Methods for achieve full tension
        %
        % 2019/02/26 - The preferred method for now is the iterated
        % interquartile Anderson-Darling test.
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
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
        
        function setSigmaFromFullTensionSolution(self)
            self.setToFullTensionWithIteratedIQAD();
            self.sigma = self.noiseDistribution.w(self.epsilon);
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
            f = @(z) abs( (alpha/(1-alpha))*outlierDistribution.pdf(abs(z))/noiseDistribution.pdf(abs(z)) - outlierOdds);
            z_outlier = fminsearch(f,sqrt(noiseDistribution.variance));
            if isempty(z_outlier)
               error('Unable to locate the outlier/noise pdf ratio requested.'); 
            end
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
        

        
        function [z,expectedVariance] = locationOfNoiseLikelihood(alpha,outlierDistribution,noiseDistribution,likelihood)
           % Finds the point above which the fewer than likelihood
           % percentage of points will be part of the noise distribution.
           %
           % For example, likelihood = 0.01 means only 1 percent of points
           % greater than z are expected to be part of the noise
           % distribution--99% are expected to be outliers.
           added = AddedDistribution(alpha,outlierDistribution,noiseDistribution);
           f = @(z) abs( ((1-alpha)*(1-noiseDistribution.cdf(z)))./(1-added.cdf(z)) - likelihood );
           z = abs(fminsearch(f,sqrt(noiseDistribution.variance)));
           expectedVariance = added.varianceInRange(-z,z);
        end
        

        
    end
end

