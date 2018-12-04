classdef GPSTensionSpline < handle
    % GPSTensionSpline Fit noisy GPS data with a tensioned interpolating
    % spline
    %   3 argument initialization
    %       f = GPSTensionSpline(t,x,y);
    %   where
    %       t       time (seconds or datetime)
    %       x       projected x-coordinate
    %       y       projected y-coordinate
    %       f       spline interpolant
    %
    %   GPSTensionSpline takes a number of optional input argument pairs.
    
    
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
        
        indicesOfOutliers
        
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
        function self = GPSTensionSpline(t,x,y,varargin)
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
            
            if shouldIdentifyOutliers == 0
                % if we're not identifying outliers, then just go ahead and
                % find the optimal expected mean square error
                self.spline_x = TensionSpline(self.t, self.x, sqrt(variance_of_the_noise), 'K', self.K, 'T', self.T, 'weightFunction', w, 'lambda', Lambda.optimalIterated);
                self.spline_y = TensionSpline(self.t, self.y, sqrt(variance_of_the_noise), 'K', self.K, 'T', self.T, 'weightFunction', w, 'lambda', Lambda.optimalIterated);
            else
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
                        outlierOdds = 10000;
                    end
                    outlierThreshold = 1/outlierOdds;
                    gps_cdf = @(z) abs(tcdf(z/sigma,nu) - outlierThreshold/2);
                    outlierDistance = fminsearch( gps_cdf, -50, optimset('TolX', 0.001, 'TolFun', 0.001) );
                end
                valid_range = [-1 1]*abs(outlierDistance);
                
                % start with a 'full tension' guess at the solution
                self.spline_x = TensionSpline(self.t, self.x, sqrt(variance_of_the_noise), 'K', self.K, 'T', self.T, 'weightFunction', w, 'lambda', Lambda.fullTensionExpected);
                self.spline_y = TensionSpline(self.t, self.y, sqrt(variance_of_the_noise), 'K', self.K, 'T', self.T, 'weightFunction', w, 'lambda', Lambda.fullTensionExpected);
                
                % now minimize the KS error with the expected distribution
                % within the valid range
                self.spline_x.Minimize( @(spline) KolmogorovSmirnovErrorForTDistribution(spline.epsilon,sigma,nu,valid_range))
                self.spline_y.Minimize( @(spline) KolmogorovSmirnovErrorForTDistribution(spline.epsilon,sigma,nu,valid_range))
                
                
            end
            
            
        end
  
    end
    
    methods (Static)
        
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
            
            gps_cdf = @(z) GPSTensionSpline.tcdf(z/sigma,nu);
            
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
