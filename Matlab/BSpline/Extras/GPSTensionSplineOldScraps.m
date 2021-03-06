classdef GPSSmoothingSpline < handle
    % GPSSmoothingSpline Fit noisy GPS data with a tensioned interpolating
    % spline
    %   3 argument initialization
    %       f = GPSSmoothingSpline(t,x,y);
    %   where
    %       t       time (seconds or datetime)
    %       x       projected x-coordinate
    %       y       projected y-coordinate
    %       f       spline interpolant
    %
    %   GPSSmoothingSpline takes a number of optional input argument pairs.
    
    
    properties
        K           % order of spline (degree = order + 1)
        T           % degree at which tension is applied
        
        % bivariate gps data
        x
        y
        t
        
        % splines
        spline_x
        spline_y
        
        distanceError = []
        indicesOfOutliers = []
        goodIndices = []
        
        % t-distribution parameters
        % variance_of_the_noise = sigma*sigma*nu/(nu-2)
        sigma = 8.5
        nu = 4.5
    end
    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Initialization
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function self = GPSSmoothingSpline(t,x,y,varargin)
            N = length(t);
            self.t = reshape(t,[],1);
            self.x = reshape(x,[],1);
            self.y = reshape(y,[],1);
            
            if length(x) ~= N || length(y) ~= N
                error('x, y, and t must have the same length.');
            end
            
            self.K = 4; % default spline order (cubic spline)
            self.T = []; % default tension *degree* (order-1)
            shouldIdentifyOutliers = 1;
            outlierOdds = []; % set odds to 1 in 10,000
            outlierDistance = [];
            
            nargs = length(varargin);
            if mod(nargs,2) ~= 0
                error('Arguments must be given as name/value pairs');
            end
                        
            for k = 1:2:length(varargin)
                if strcmp(varargin{k}, 'K')
                    self.K = varargin{k+1};
                elseif strcmp(varargin{k}, 'S')
                    self.K = varargin{k+1}+1;
                elseif strcmp(varargin{k}, 'T')
                    self.T = varargin{k+1};
                elseif strcmp(varargin{k}, 'sigma')
                    self.sigma = varargin{k+1};
                elseif strcmp(varargin{k}, 'nu')
                    self.nu = varargin{k+1};
                elseif strcmp(varargin{k}, 'shouldIdentifyOutliers')
                    shouldIdentifyOutliers = varargin{k+1};
                elseif strcmp(varargin{k}, 'outlierOdds')
                    outlierOdds = varargin{k+1};
                elseif strcmp(varargin{k}, 'outlierDistance')
                    outlierDistance = varargin{k+1};
                end
            end
            
            
            % initialization complete, let's do some work
            sigma = self.sigma; nu = self.nu;
            variance_of_the_noise = sigma*sigma*nu/(nu-2);
            w = @(z)((nu/(nu+1))*sigma^2*(1+z.^2/(nu*sigma^2)));
            
            if shouldIdentifyOutliers == 1
                % we're looking for outliers, that is, data points that
                % have low odds of occurring. By default we will find
                % outliers with less than 1 in 10,000 odds of occurring,
                % and then convert that to a distance.
                %
                % Note that you can just match the sample variance with the
                % expected variance of the noise because outliers (by
                % definition!) will add too much variance.
                if isempty(outlierDistance)
                    if isempty(outlierOdds)
                        outlierOdds = length(t);
                    end
                    outlierThreshold = 0.001;
                    gps_cdf = @(z) abs(tcdf(z/sigma,nu) - outlierThreshold/2);
                    outlierDistance = fminsearch( gps_cdf, -50, optimset('TolX', 0.001, 'TolFun', 0.001) );                   
%                     gps_cdf = @(z) abs(tcdf(z/sigma,nu).^length(t) - .75);
%                     outlierDistance = fminsearch( gps_cdf, 50, optimset('TolX', 0.001, 'TolFun', 0.001) );

%                     gps_cdf = @(z) abs( (2*tcdf(abs(z/sigma),nu)-1).^length(t) - 0.95);
%                     outlierDistance = fminsearch( gps_cdf, 50, optimset('TolX', 0.001, 'TolFun', 0.001) );
                end
                valid_range = [-1 1]*abs(outlierDistance);
                
                % start with a 'full tension' guess at the solution
                self.spline_x = SmoothingSpline(self.t, self.x, sqrt(variance_of_the_noise), 'K', self.K, 'T', self.T, 'weightFunction', w, 'lambda', Lambda.fullTensionExpected);
                self.spline_y = SmoothingSpline(self.t, self.y, sqrt(variance_of_the_noise), 'K', self.K, 'T', self.T, 'weightFunction', w, 'lambda', Lambda.fullTensionExpected);
                
                % now minimize the KS error with the expected distribution
                % within the valid range (i.e., excluding outliers)
                fprintf('Finding full tension solution by ignoring points with errors more than %.1f meters.\n',abs(outlierDistance));
                self.spline_x.Minimize( @(spline) GPSSmoothingSpline.KolmogorovSmirnovErrorForTDistribution(spline.epsilon,sigma,nu,valid_range));
                self.spline_y.Minimize( @(spline) GPSSmoothingSpline.KolmogorovSmirnovErrorForTDistribution(spline.epsilon,sigma,nu,valid_range));
                
                self.distanceError = sqrt( self.spline_x.epsilon.^2 + self.spline_y.epsilon.^2 );
                
                [r, ~, cdf] = GPSSmoothingSpline.TwoDimStudentTProbabilityDistributionFunction(sigma, nu);
%                 outlierCut = interp1(cdf,r,1-outlierThreshold,'spline'); % set to spline so that it extrapolates by default, which is better than returning NaN.
                
                % About 10% of the time we will have one legit point above
                % this threshold.
                outlierThreshold = 0.05; % p-value test---control false positives
                jpd = cdf.^length(t);
                idx = jpd>0.1; % make things monotonic (we don't need that data anyway)
                outlierCut = interp1(jpd(idx),r(idx),1-outlierThreshold,'spline');
                
                self.indicesOfOutliers = find(self.distanceError(2:end-1) >= outlierCut)+1;
                
                fprintf('Found %d points (of %d) that exceed the two-dimensional outlier distance of %.1f meters.\n',length(self.indicesOfOutliers), length(self.t), outlierCut); 
                fprintf('This threshold was chosen because you would expect only 1 point to exceed this distance %d%% of the time.\n',round(outlierThreshold*100));
            end
            
            self.goodIndices = setdiff((1:length(t))',self.indicesOfOutliers);
            
            self.spline_x = SmoothingSpline(self.t(self.goodIndices), self.x(self.goodIndices), sqrt(variance_of_the_noise), 'K', self.K, 'T', self.T, 'weightFunction', w, 'lambda', Lambda.optimalIterated);
            self.spline_y = SmoothingSpline(self.t(self.goodIndices), self.y(self.goodIndices), sqrt(variance_of_the_noise), 'K', self.K, 'T', self.T, 'weightFunction', w, 'lambda', Lambda.optimalIterated);         
        end
        
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
                    
                    varargout{1} = self.spline_x(time,NumDerivatives);
                    varargout{2} = self.spline_y(time,NumDerivatives);
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GET %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                case '.'
                    [varargout{1:nargout}] = builtin('subsref',self,index);
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RESTRICT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                case '{}'
                    error('The GPSSmoothingSpline class does not know what to do with {}.');
                otherwise
                    error('Unexpected syntax');
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

        function lambda = MinimizeFunctionOfTwoSplines(aSmoothingSpline1,aSmoothingSpline2,functionOfSplines)
            epsilon = 1e-15;
            errorFunction = @(log10lambdaPlusEpsilon) GPSSmoothingSpline.FunctionOfSplineWrapper(aSmoothingSpline1,aSmoothingSpline2,log10lambdaPlusEpsilon,functionOfSplines);
            optimalLog10lambdaPlusEpsilon = fminsearch( errorFunction, log10(aSmoothingSpline1.lambda+epsilon), optimset('TolX', 0.01, 'TolFun', 0.01) );
            lambda = 10^optimalLog10lambdaPlusEpsilon - epsilon;
        end
        
        function error = FunctionOfSplineWrapper(aSmoothingSpline1, aSmoothingSpline2, log10lambdaPlusEpsilon, functionOfSplines)
            epsilon = 1e-15;
            aSmoothingSpline1.lambda = 10^log10lambdaPlusEpsilon-epsilon;  
            aSmoothingSpline2.lambda = 10^log10lambdaPlusEpsilon-epsilon;  
            error = functionOfSplines(aSmoothingSpline1,aSmoothingSpline2);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % 2D LogLikelihood for Student-t and Normal
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function likelihood = LogLikelihood2D(spline_x, spline_y, sigma_t, nu, sigma_g)
                        
%             distance_pdf = @(epsilon) interp1(r,pdf,epsilon,'spline');
%             acceleration_pdf = @(r) r.*exp(-r.*r/(2*sigma_T*sigma_T))/(sigma_T*sigma_T);

            [r, pdf] = GPSSmoothingSpline.TwoDimStudentTProbabilityDistributionFunction(sigma_t, nu);
            pdf(pdf<1e-15)=1e-15; % avoid log(0).
            log_distance_pdf = @(epsilon) interp1(r,log(pdf),epsilon,'spline');
            log_acceleration_pdf = @(r) log(r) + (-r.*r/(2*sigma_g*sigma_g));
            
            distance = sqrt(spline_x.epsilon.^2 + spline_y.epsilon.^2);
            acceleration = sqrt( spline_x.UniqueValuesAtHighestDerivative().^2 + spline_y.UniqueValuesAtHighestDerivative().^2);
            
            error_likelihood = -sum(log_distance_pdf(distance));
            acceleration_likelihood = -sum(log_acceleration_pdf(acceleration));
            
            likelihood = error_likelihood + acceleration_likelihood;
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % 
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function totalError = KolmogorovSmirnovErrorForTDistribution( epsilon, sigma, nu, range)
            %%KolmogorovSmirnovErrorForTDistribution
            %
            % Return a measure of the Kolmogorov-Smirnov test metric for the
            % distribution matched to the position errors
            %
            % the range allows you to restrict the range (in meters) over which the
            % test is applied.
            
            gps_cdf = @(z) GPSSmoothingSpline.tcdf(z/sigma,nu);
            
            if nargin == 4
                x = sort(epsilon( epsilon > range(1) & epsilon < range(2) ));
            else
                x = sort(epsilon);
            end
            
            n = length(x);
            y_data = (1:n)'/n;
            y = gps_cdf(x);
            
            D = max(abs(y-y_data));
            
            totalError = (sqrt(n) + 0.12 + 0.11/sqrt(n))*D;
        end
        
        function [r, pdf, cdf] = TwoDimStudentTProbabilityDistributionFunction(sigma, nu)
            % returns the pdf and cdf of the two-dimensional t-distribution
            % as a function of distance, r.
            
            % This function needs to fast, and the cdf needs to be
            % accurate. So we use a strategy of putting points where
            % necessary (logarthimically placed), and using the cdf as the
            % fundamental quantity (rather than the pdf) because we can
            % adjust for discretization more easily (b/c it's already
            % integrated).
            
            % scale our axis to be 20 times the standard deviation of the
            % 1 dimensional pdf.
            maxR = 20*sqrt(sigma*sigma*nu/(nu-2));
            
            % This appears to be big enough, although we can go bigger if
            % needed.
            N = 150;
            
            % create a logarithmic axis that includes zero
            r = [0;10.^(linspace(log10(maxR/1e3),log10(maxR),N))'];
            x = [-flip(r(2:end),1); r];
            
            % evaluate the cdf on that axis
            pdf_int = tcdf(x/sigma,nu);
            
            % use that to define the pdf...
            dx = diff(x);
            pdf = diff(pdf_int)./dx;
            %...so that sum(pdf.*dx)=1.00000
            
            % now create a 2d version
            pdf2d_norm = (pdf .* reshape(pdf,1,[])) .* (dx .* reshape(dx,1,[]));
            % again, where sum(pdf2d_norm(:)) = 1.0000
            
            % create a 2d distance metric, rho
            y = [-flip(r(2:end)); r(2:end)];
            rho = sqrt( y.*y + reshape(y.*y,1,[]));
            
            % we are going to bin using the diagonal of rho, call it s.
            s = diag(rho);
            s = [0; s((N+1):end)];
            
            % Now sum up the energy in each bin
            pdf1d = zeros( size(s) );
            midS = [0; s(1:(end-1))+diff(s)/2];
            for i = 1:(length(s)-1)
                pdf1d(i) = sum( pdf2d_norm(  rho >= midS(i) & rho <= midS(i+1) ) );
            end
            % it must be true that sum(pdf1d) = 1.0000
            
            % but we want that sum(pdf1d .* diff(s)) = 1.000
            pdf1d(2:end) = pdf1d(2:end)./diff(s);
            cdf1d = [0; cumsum(pdf1d(2:end).*diff(s))];
            
            % now take care of the the fact that we used the diagnoal of a
            % 2d rectangle, rather than the inscribed circle.
            
            pdf = pdf1d(s<maxR);
            cdf = cdf1d(s<maxR);
            r = s(s<maxR);
        end
        
        function p = tcdf(x,n)
            % TCDF returns student cumulative distribtion function
            %
            % cdf = tcdf(x,DF);
            %
            % Computes the CDF of the students distribution
            %    with DF degrees of freedom
            % x,DF must be matrices of same size, or any one can be a scalar.
            %
            % see also: NORMCDF, NORMPDF, NORMINV
            
            % Reference(s):
            
            %	$Revision: 1.1 $
            %	$Id: tcdf.m,v 1.1 2003/09/12 12:14:45 schloegl Exp $
            %	Copyright (c) 2000-2003 by  Alois Schloegl <a.schloegl@ieee.org>
            
            %    This program is free software; you can redistribute it and/or modify
            %    it under the terms of the GNU General Public License as published by
            %    the Free Software Foundation; either version 2 of the License, or
            %    (at your option) any later version.
            %
            %    This program is distributed in the hope that it will be useful,
            %    but WITHOUT ANY WARRANTY; without even the implied warranty of
            %    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
            %    GNU General Public License for more details.
            %
            %    You should have received a copy of the GNU General Public License
            %    along with this program; if not, write to the Free Software
            %    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
            
            
            % check size of arguments
            n = x+n-x;	  % if this line causes an error, size of input arguments do not fit.
            z = n ./ (n + x.^2);
            
            % allocate memory
            p = z;
            p(x==Inf) = 1;
            
            % workaround for invalid arguments in BETAINC
            tmp   = isnan(z) | ~(n>0);
            p(tmp)= NaN;
            ix    = (~tmp);
            p(ix) = betainc (z(ix), n(ix)/2, 1/2) / 2;
            
            ix    = find(x>0);
            p(ix) = 1 - p(ix);
            
            % shape output
            p = reshape(p,size(z));
        end
        
    end
end
