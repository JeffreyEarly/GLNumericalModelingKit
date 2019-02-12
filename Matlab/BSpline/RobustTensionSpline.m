classdef RobustTensionSpline < TensionSpline
    %RobustTensionSpline Allows for outliers
    
    properties
        noiseDistribution
        outlierDistribution
        
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
            self.outlierThreshold = self.noiseDistribution.locationOfCDFPercentile(1-1/10000/2);

            
%             self.secondIteration();
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
        
        function setToFullTensionWithInnerSV(self,alpha)
            % This algorithm takes the interquartile range, and requires it
            % match the expected variance in that range. This is different than
            % above in that it has a fixed number of points.
            
            zmin = self.distribution.locationOfCDFPercentile(alpha/2);
            zmax = self.distribution.locationOfCDFPercentile(1-alpha/2);
            expectedVarianceInRange = self.distribution.varianceInRange(zmin,zmax);
            self.minimize(@(spline) abs(spline.varianceOfInterquartile(alpha) - expectedVarianceInRange));
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
        
        function [outlierDistribution, alpha,zOutlier] = estimateOutlierDistribution(self)
            self.setToFullTensionWithSV(0.5);
            
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
        
        function [outlierDistribution, alpha,zOutlier] = estimateOutlierDistributionABetterMethod(self)
            self.setToFullTensionWithInnerSV(0.5);
            
            % Let's find a reasonable set of z_crossover points.
            epsilon = self.epsilon;
            abs_eps = sort(abs(epsilon),'descend');
            n = length(epsilon);
            cdf_eps = (1:n)'/n;
            
            lastIndex = find( cdf_eps > 0.6, 1, 'first');
            
            if lastIndex < 10
                error('Did not expect this')
            end
            
            s2_total = mean(epsilon.^2);
            s2_noise = self.noiseDistribution.variance;
            
            if s2_total/s2_noise < 1.1
                outlierDistribution = [];
                alpha = 0;
                return;
            end
            
            alpha_outlier = 10.^(linspace(log10(0.01),log10(0.5),100))';
            sigma_outlier = (s2_total-(1-alpha_outlier)*s2_noise)./(3*alpha_outlier);
            
            ks_error = zeros(size(alpha_outlier));
            
            % handle the alpha=0 case separately.
            ks_error(1) = self.noiseDistribution.kolmogorovSmirnovError(epsilon);
            
            for iAlpha = 1:length(alpha_outlier)
                newAddedDistribution = AddedDistribution(alpha_outlier(iAlpha),StudentTDistribution(sqrt(sigma_outlier(iAlpha)),3),self.noiseDistribution);
                ks_error(iAlpha) = newAddedDistribution.kolmogorovSmirnovError(epsilon);   
            end
            
            [minAlpha,minIndex] = min(ks_error);
            if (self.noiseDistribution.kolmogorovSmirnovError(epsilon) < minAlpha)
                outlierDistribution = [];
                alpha = 0;
            else
                alpha = alpha_outlier(minIndex);
                sigma2o = sigma_outlier(minIndex);
                outlierDistribution = StudentTDistribution(sqrt(sigma2o),3);
                
                f = @(z) abs( (1-alpha)*self.noiseDistribution.pdf(-abs(z)) - alpha*outlierDistribution.pdf(-abs(z)) );

                zOutlier = abs(fminsearch(f,sqrt(s2_noise)));
            end
        end
        
        
        function rebuildOutlierDistribution(self)
            [newOutlierDistribution, alpha,zOutlier] = self.estimateOutlierDistributionABetterMethod();
            
            if alpha > 0.0
                newAddedDistribution = AddedDistribution(alpha,newOutlierDistribution,self.noiseDistribution);
                self.distribution = newAddedDistribution;
                [newOutlierDistribution, alpha,zOutlier] = self.estimateOutlierDistributionABetterMethod();
                
                crossover = @(z) abs((1-alpha).*self.noiseDistribution.pdf(z) - alpha.*newOutlierDistribution.pdf(z));
                z_crossover = abs(fminsearch(crossover,-abs(zOutlier)));
                
                epsilon = self.epsilon;
                noiseIndices = epsilon >= -z_crossover & epsilon <= z_crossover;
                outlierIndices = ~noiseIndices;
                
                self.distribution = newAddedDistribution;
                self.distribution.w = @(z) noiseIndices .* self.noiseDistribution.w(z) + outlierIndices .* newOutlierDistribution.w(z);
                
                self.minimize( @(spline) spline.expectedMeanSquareErrorInRange(-abs(z_crossover),z_crossover) );
            else
                self.minimizeExpectedMeanSquareError();
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

