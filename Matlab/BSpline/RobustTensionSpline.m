classdef RobustTensionSpline < TensionSpline
    %RobustTensionSpline Allows for outliers
    
    properties
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
            % only 1 percent outliers.
            
%             nu = 3.0;
%             noiseDistribution = StudentTDistribution(sqrt(distribution.variance*(nu-2)/nu),nu);
%             
%             distribution =  noiseDistribution;
            

            
            % Override any user settings for lambda
            didOverrideLambda = 0;
            alpha = 0.01;
            for k = 1:2:length(varargin)
                if strcmp(varargin{k}, 'lambda')
                    varargin{k+1} = Lambda.fullTensionExpected;
                    didOverrideLambda = 1;
                elseif strcmp(varargin{k}, 'alpha')
                    alpha = varargin{k+1};
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
        end
        
        function firstIteration(self,alpha)
            % Minimize using the expected mean square error
            self.zmin = self.noiseDistribution.locationOfCDFPercentile(alpha/2);
            self.zmax = self.noiseDistribution.locationOfCDFPercentile(1-alpha/2);
            self.minimize( @(spline) spline.expectedMeanSquareErrorInRange(self.zmin,self.zmax) );
            
            % Remove knot support from outliers, and rescale the
            % distribution
            self.indicesOfOutliers = find(abs(self.epsilon) > self.outlierThreshold);
        end
        
        function firstIterationCV(self)
            % Minimize using the expected mean square error from
            % cross-validation
            self.minimize( @(spline) spline.expectedMeanSquareErrorFromCV() );
            
            self.indicesOfOutliers = find(abs(self.epsilon) > self.outlierThreshold);
        end
        
        function setToFullTensionWithSV(self,alpha)
            % This gives us "full tension" using the sample variance.
            if nargin == 1
                alpha = 1/10;
            end
%             self.minimize(@(spline) abs(spline.sampleVarianceInPercentileRange(alpha/2,1-alpha/2)-spline.noiseDistribution.varianceInPercentileRange(alpha/2,1-alpha/2)));
%             self.minimize(@(spline) abs(spline.sampleVarianceInPercentileRange(alpha/2,1-alpha/2)-(1-1/spline.effectiveSampleSizeFromVarianceOfTheMeanInPercentileRange(alpha/2,1-alpha/2))*spline.noiseDistribution.varianceInPercentileRange(alpha/2,1-alpha/2)));

            zmin = self.noiseDistribution.locationOfCDFPercentile(alpha/2);
            zmax = self.noiseDistribution.locationOfCDFPercentile(1-alpha/2);
            expectedVarianceInRange = self.noiseDistribution.varianceInRange(zmin,zmax);
            self.minimize(@(spline) abs(spline.sampleVarianceInRange(zmin,zmax) - expectedVarianceInRange));
        end
        
        function setToFullTensionWithInnerSVOnNoiseDistribution(self,alpha)
            % This algorithm takes the interquartile range, and requires it
            % match the expected variance in that range. This is different than
            % above in that it has a fixed number of points.
            
            zmin = self.noiseDistribution.locationOfCDFPercentile(alpha/2);
            zmax = self.noiseDistribution.locationOfCDFPercentile(1-alpha/2);
            expectedVarianceInRange = self.noiseDistribution.varianceInRange(zmin,zmax);
            self.minimize(@(spline) abs(spline.varianceOfInterquartile(alpha) - expectedVarianceInRange));
        end
        
        function setToFullTensionWithInnerSV(self,alpha)
            % This algorithm takes the interquartile range, and requires it
            % match the expected variance in that range. This is different than
            % above in that it has a fixed number of points.
            
            zmin = self.distribution.locationOfCDFPercentile(alpha/2);
            zmax = self.distribution.locationOfCDFPercentile(1-alpha/2);
            expectedVarianceInRange = self.distribution.varianceInRange(zmin,zmax);
            self.setInnerVarianceToExpectedValue(alpha,expectedVarianceInRange);
        end
        
        function setInnerVarianceToExpectedValue(self,alpha,expectedVarianceInRange)
            fprintf('Setting variance to %.1f m^2\n',expectedVarianceInRange);
            self.minimize(@(spline) abs(spline.varianceOfInterquartile(alpha) - expectedVarianceInRange));
        end
        
        function setToFullTensionWithIteratedIQSV(self,extraVarianceFactor)
           % set to full tension by matching the expected sample variance with the interquartile sample variance.
           % at each iteration the outlier distribution is estimated
           if nargin < 2
               extraVarianceFactor = 1.0;
           end
           pctRange = 1/2; % alpha=1/2 is the inner half of the distribution
           self.setInnerVarianceToExpectedValue(pctRange,extraVarianceFactor*self.noiseDistribution.varianceInPercentileRange(pctRange/2,1-pctRange/2));
           lastAlpha = 0.0;
           totalIterations = 0;
           [newOutlierDistribution, newAlpha] = self.estimateOutlierDistribution();
           while (abs(lastAlpha-newAlpha) > 0.01 && totalIterations < 10)
               addedDistribution = AddedDistribution(newAlpha,newOutlierDistribution,self.noiseDistribution);
               self.setInnerVarianceToExpectedValue(pctRange,extraVarianceFactor*addedDistribution.varianceInPercentileRange(pctRange/2,1-pctRange/2));
               lastAlpha = newAlpha;
               [newOutlierDistribution, newAlpha] = self.estimateOutlierDistribution();
               totalIterations = totalIterations + 1;
           end
        end
        
        function setToFullTensionWithIteratedIQAD(self)
            % set to full tension by minimizing the Anderson-Darling test
            % on the interquartile range of epsilons.
            self.minimize(@(spline) self.noiseDistribution.andersonDarlingInterquartileError(spline.epsilon));
            lastAlpha = 0.0;
            totalIterations = 0;
            [newOutlierDistribution, newAlpha] = self.estimateOutlierDistribution();
            while (abs(lastAlpha-newAlpha) > 0.01 && totalIterations < 10)
                if newAlpha > 0 && ~isempty(newOutlierDistribution)
                    addedDistribution = AddedDistribution(newAlpha,newOutlierDistribution,self.noiseDistribution);
                else
                    addedDistribution = self.noiseDistribution;
                end
                self.minimize(@(spline) addedDistribution.andersonDarlingInterquartileError(spline.epsilon));
                lastAlpha = newAlpha;
                [newOutlierDistribution, newAlpha] = self.estimateOutlierDistribution();
                totalIterations = totalIterations + 1;
            end
        end
        
        function generateEpsilonWeighting(self)
           % assuming we're at full tension
           if self.alpha > 0 && ~isempty(self.outlierDistribution)
               epsilon = self.epsilon;
               self.w_epsilon = (1-self.alpha)*self.noiseDistribution.pdf(epsilon)./self.distribution.pdf(epsilon);
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
        
        
        function sampleVariance = varianceOfInterquartile(self,alpha)
           epsilon = sort(self.epsilon);
           indices = floor(length(epsilon)*alpha/2):ceil(length(epsilon)*(1-alpha/2));
           sampleVariance = mean(epsilon(indices).^2);
        end
        
        function setToFullTensionWithKS(self,alpha)
            % This gives us "full tension" using a Kolmogorov-Smirnov test.
            if nargin == 1
                alpha = 1/10;
            end
            
            epsilon_min = self.noiseDistribution.locationOfCDFPercentile(alpha/2);
            epsilon_max = self.noiseDistribution.locationOfCDFPercentile(1-alpha/2);            
            self.minimize( @(spline) self.noiseDistribution.kolmogorovSmirnovError(spline.epsilon,epsilon_min,epsilon_max) );
        end
        
        function [outlierDistribution, alpha,zOutlier] = estimateOutlierDistributionOldMethod(self)
%             self.setToFullTensionWithSV(0.5);
            
            % Identify the spot where our distribution no longer looks
            % right using the KS test
            epsilon = self.epsilon;
            z = linspace(self.noiseDistribution.locationOfCDFPercentile(1-1/100000/2),self.noiseDistribution.locationOfCDFPercentile(1-1/5/2),10)';
            err = zeros(length(z),1);
            for iZ = 1:length(z)
                err(iZ) = self.noiseDistribution.kolmogorovSmirnovError(epsilon,-abs(z(iZ)),abs(z(iZ)));
            end
            [~,iZ] = min(err);
            minZ = z(iZ);
            zOutlier = fminsearch( @(z) self.noiseDistribution.kolmogorovSmirnovError(epsilon,-abs(z),abs(z)), minZ);
%             fprintf('zOutlier: %f\n',zOutlier);
            
            % Now try to estimate the outlier distribution
            noiseIndices = find(epsilon >= -zOutlier & epsilon <= zOutlier);
            outlierIndices = setdiff(1:length(self.t),noiseIndices);
            
            n_eff_o = self.effectiveSampleSizeFromVarianceOfTheMeanForIndices(outlierIndices);
            n_eff_n = self.effectiveSampleSizeFromVarianceOfTheMeanForIndices(noiseIndices);
            s2_total = mean(epsilon.^2);
            s2_noise = (1-1/n_eff_n)*self.noiseDistribution.variance;
            
            if s2_total/s2_noise < 1.2
                outlierDistribution = [];
                alpha = 0;
                return;
            end
            
            % if s2_total is anywhere near the expected variance, we should
            % bail
            
            alpha = @(sigma2) (s2_total-s2_noise)/(3*(1-1/n_eff_o)*sigma2 - s2_noise);
            
            outlierFraction = length(outlierIndices)/length(self.t);
            
            f = @(sigma2) abs( (1-alpha(abs(sigma2)))*2*self.noiseDistribution.cdf(-abs(zOutlier)) + alpha(abs(sigma2))*2*StudentTDistribution(sqrt(abs(sigma2)),3).cdf(-abs(zOutlier)) - outlierFraction);
            
            sigma2_outlier = fminsearch(f,s2_total);
            alpha = alpha(sigma2_outlier);
            
            % if sigma2_outlier is less than sigma2_noise, we should bail.
            
            outlierDistribution = StudentTDistribution(sqrt(sigma2_outlier),3);
        end
        
        function [outlierDistribution, alpha,zOutlier] = estimateOutlierDistribution(self)
%             self.setToFullTensionWithInnerSV(0.5);
            
            % Let's find a reasonable set of z_crossover points.
            epsilon = self.epsilon;            
            s2_total = mean(epsilon.^2);
            s2_noise = self.noiseDistribution.variance;
            
            % What's a reasonable cutoff here?
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
                newAddedDistribution = AddedDistribution(alpha_outlier(iAlpha),StudentTDistribution(sqrt(sigma2_outlier(iAlpha)),3),self.noiseDistribution);
                var_total(iAlpha) = newAddedDistribution.variance;
                ks_error(iAlpha) = newAddedDistribution.andersonDarlingError(epsilon);   
            end
            
            [minError,minIndex] = min(ks_error);
            if (self.noiseDistribution.andersonDarlingError(epsilon) < minError)
                outlierDistribution = [];
                alpha = 0;
                zOutlier = Inf;
            else
                alpha = alpha_outlier(minIndex);
                sigma2o = sigma2_outlier(minIndex);
                outlierDistribution = StudentTDistribution(sqrt(sigma2o),3);
                
                if nargout == 3
                    f = @(z) abs( (1-alpha)*self.noiseDistribution.pdf(z) - alpha*outlierDistribution.pdf(z) );
                    zOutlier = abs(fminsearch(f,sqrt(s2_noise)));
                end
            end
        end
        
        
        function rebuildOutlierDistribution(self)
            [self.outlierDistribution, self.alpha] = self.estimateOutlierDistribution();
            
            if self.alpha > 0.0
                fprintf('Rebuilding outlier distribution with alpha=%.2f and sqrt(var)=%.1f\n',self.alpha,sqrt(self.outlierDistribution.variance));
                self.distribution = AddedDistribution(self.alpha,self.outlierDistribution,self.noiseDistribution);
            end
        end
        
        function z_crossover = rebuildOutlierDistributionAndAdjustWeightings(self,outlierOdds)
            [self.outlierDistribution, self.alpha,z_crossover] = self.estimateOutlierDistribution();
            
            if self.alpha > 0.0
                fprintf('Rebuilding outlier distribution/weightings with alpha=%.2f and sqrt(var)=%.1f\n',self.alpha,sqrt(self.outlierDistribution.variance));
                self.distribution = AddedDistribution(self.alpha,self.outlierDistribution,self.noiseDistribution);
                
                if ~isempty(outlierOdds)
                   f = @(z) abs( (self.alpha/(1-self.alpha))*self.outlierDistribution.pdf(-abs(z))/self.noiseDistribution.pdf(-abs(z)) - outlierOdds);
                   z_outlier = fminsearch(f,sqrt(self.noiseDistribution.variance));
                   fprintf('Setting outlier cutoff at z=%.1f m\n',z_outlier);
                else
                    z_outlier = z_crossover;
                end
                
                noiseIndices = abs(self.epsilon) <= z_outlier;
                
%                 self.distribution.w = @(z) noiseIndices .* self.noiseDistribution.w(z) + (~noiseIndices) .* self.outlierDistribution.w(z);
                self.sigma = noiseIndices .* sqrt(self.noiseDistribution.variance) + (~noiseIndices) .* sqrt(self.outlierDistribution.variance);

                self.distribution = self.noiseDistribution;
%                 self.minimize( @(spline) spline.expectedMeanSquareErrorInRange(-abs(z_crossover),z_crossover) );
%             else
%                 self.minimizeExpectedMeanSquareError();
            end
%             fprintf('z_crossover: %.2f\n',z_crossover);
%             self.removeOutlierKnotsAndRetensionInRange(-abs(z_crossover),z_crossover);
        end
        
        function rebuildOutlierDistributionAndAdjustWeightingsOldMethod(self)
            [newOutlierDistribution, alpha,z_crossover] = self.estimateOutlierDistributionOldMethod();
            
            if alpha > 0.0
                fprintf('Rebuilding outlier distribution/weightings with alpha=%.2f and sqrt(var)=%.1f\n',alpha,sqrt(newOutlierDistribution.variance));
                newAddedDistribution = AddedDistribution(alpha,newOutlierDistribution,self.noiseDistribution);
                self.distribution = newAddedDistribution;
                
                epsilon = self.epsilon;
                noiseIndices = epsilon >= -z_crossover & epsilon <= z_crossover;
                outlierIndices = ~noiseIndices;
                
                self.distribution = newAddedDistribution;
                self.distribution.w = @(z) noiseIndices .* self.noiseDistribution.w(z) + outlierIndices .* newOutlierDistribution.w(z);
                
                %                 self.minimize( @(spline) spline.expectedMeanSquareErrorInRange(-abs(z_crossover),z_crossover) );
                %             else
                %                 self.minimizeExpectedMeanSquareError();
            end
            %             fprintf('z_crossover: %.2f\n',z_crossover);
            %             self.removeOutlierKnotsAndRetensionInRange(-abs(z_crossover),z_crossover);
        end

        function rescaleDistributionAndRetension(self,alpha)
            % you could try to scale the variance correctly.
            % we have a certain percentage of points outside zmin/zmax
            % we have a certain variance of those points
            % is that enough to properly rescale an added distribution?
            % I want a dist w/ x^2 variance in a certain range
            % I then want the added cdf to give the right sum at zmin/zmax
            
            % Remove knot support from outliers, and rescale the
            % distribution
            if isempty(self.indicesOfOutliers)
                % No outliers? Then revert to the usual case
                self.distribution = self.noiseDistribution;
            else
                % otherwise rescale the distribution more appropriately
                nu = 3.0;
                epsilon = self.epsilon;
                variance = var(epsilon(self.indicesOfOutliers));
                sigma = sqrt(variance*(nu-2)/nu);
                self.outlierDistribution = StudentTDistribution(sigma,nu);
                self.distribution = AddedDistribution(length(self.indicesOfOutliers)/length(self.t),self.outlierDistribution,self.noiseDistribution);
            end
                        
            % Minimize using the expected mean square error
            self.zmin = self.noiseDistribution.locationOfCDFPercentile(alpha/2);
            self.zmax = self.noiseDistribution.locationOfCDFPercentile(1-alpha/2);
            self.minimize( @(spline) spline.expectedMeanSquareErrorInRange(self.zmin,self.zmax) );
            self.indicesOfOutliers = find(abs(self.epsilon) > self.outlierThreshold);
        end
        
        function removeOutlierKnotsAndRetension(self,alpha)  
            newKnotIndices = self.goodIndices;
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
            
            % Minimize using the expected mean square error
            self.zmin = self.noiseDistribution.locationOfCDFPercentile(alpha/2);
            self.zmax = self.noiseDistribution.locationOfCDFPercentile(1-alpha/2);
            self.minimize( @(spline) spline.expectedMeanSquareErrorInRange(self.zmin,self.zmax) );
            self.indicesOfOutliers = find(abs(self.epsilon) > self.outlierThreshold);
        end
        
        function removeOutlierKnotsAndRetensionInRange(self,zmin,zmax)
            newKnotIndices = self.goodIndices;
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
            
            % Minimize using the expected mean square error
            self.zmin = zmin;
            self.zmax = zmax;
            self.minimize( @(spline) spline.expectedMeanSquareErrorInRange(self.zmin,self.zmax) );
            self.indicesOfOutliers = find(abs(self.epsilon) > self.outlierThreshold);
        end
        
        function minimizeExpectedMeanSquareError(self)
            % This value of beta is chosen from
            % MakeTableOutliersWithAddedDistributionBlind.m as a good value
            % for distribution with very few outliers.
            beta = 1/400;
            zmin_ = self.noiseDistribution.locationOfCDFPercentile(beta/2);
            zmax_ = self.noiseDistribution.locationOfCDFPercentile(1-beta/2);
            
            if ~(isempty(self.outlierDistribution) || self.alpha < 0.01)
                % If we a non-trivial outlier distribution, includepoints
                % that are more likely than not to be part of the
                % characterized noise distribution.
                noiseOdds = 3; % This is from MakeTableOutliersBlindMinimization.m
                f = @(z) abs( (1-self.alpha)*self.noiseDistribution.pdf(z) - noiseOdds*self.alpha*self.outlierDistribution.pdf(z) );
                zoutlier = abs(fminsearch(f,sqrt(self.noiseDistribution.variance)));     
                if zoutlier < zmax_
                    zmin_ = -zoutlier;
                    zmax_ = zoutlier;
                end
            end
            self.minimize( @(spline) spline.expectedMeanSquareErrorInRange(zmin_,zmax_) );
%             fprintf('Minimizing in range: (%.1f, %.1f) with expected variance %.1f m^2\n',zmin_,zmax_,self.distribution.varianceInRange(zmin_,zmax_));
        end
        
    end
    
    methods (Static)
        
        function [outlierDistribution, alpha] = estimateOutlierDistributionFromKnownNoise(epsilon,noiseDistribution)
            
            % Let's find a reasonable set of z_crossover points.
            s2_total = mean(epsilon.^2);
            s2_noise = noiseDistribution.variance;
            
            if s2_total/s2_noise < 1.1
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
        
        function w_epsilon = generateEpsilonWeightingFromOutlierDistribution(epsilon,alpha,outlierDistribution,noiseDistribution)
            if alpha > 0 && ~isempty(outlierDistribution)
                w_epsilon = (1-alpha)*noiseDistribution.pdf(epsilon)./AddedDistribution(alpha,outlierDistribution,noiseDistribution).pdf(epsilon);
            else
                w_epsilon = ones(size(epsilon));
            end
        end
        
    end
end

