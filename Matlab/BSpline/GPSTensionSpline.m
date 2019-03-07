classdef GPSTensionSpline < handle
    % GPSTensionSpline Fit noisy GPS data with a tensioned interpolating
    % spline
    %   3 argument initialization
    %       f = GPSTensionSpline(t,lat,lon);
    %   where
    %       t       time (seconds or datetime)
    %       lat     latitude
    %       lon     longitude
    %       f       spline interpolant
    %
    %   GPSTensionSpline takes a number of optional input argument pairs.
    
    
    properties
        K           % order of spline (degree = order + 1)
        T           % degree at which tension is applied
        
        % bivariate gps data
        t_raw % format that was given
        t % converted to second from start date
        lat
        lon
        
        % projected data
        lon0 % reference longitude of transverse mercator projection being used
        x0  % false easting
        y0  % false northing
        x
        y
        
        % translated and rotated data
        xbar
        ybar
        theta
        q
        r
        
        shouldUseRobustFit = 0;
        
        distribution
        
        % lookup table for tension
        tensionLambdaTable = [] % size(n,3) [x_T lambda_x lambda_y];
        
        % splines
        spline_x
        spline_y
        
        % t-distribution parameters
        % variance_of_the_noise = sigma*sigma*nu/(nu-2)
        sigma = 8.5
        nu = 4.5
        
        k0 = 0.9996
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
        function self = GPSTensionSpline(t,lat,lon,varargin)
            N = length(t);
            self.t_raw = reshape(t,[],1);
            self.lat = reshape(lat,[],1);
            self.lon = reshape(lon,[],1);
            
            if isdatetime(self.t_raw)
                self.t = (datenum(self.t_raw) - datenum(self.t_raw(1)))*24*60*60;
            else
                self.t = self.t_raw;
            end
            
            if length(lat) ~= N || length(lon) ~= N
                error('lat, lon, and t must have the same length.');
            end
            
            self.K = 4; % default spline order (cubic spline)
            self.T = self.K-1; % default tension *degree* (order-1)
            self.lon0 = min(lon)+(max(lon)-min(lon))/2;
            self.x0 = NaN;
            self.y0 = NaN;
            
            
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
                elseif strcmp(varargin{k}, 'lon0')
                    self.lon0 = varargin{k+1};
                elseif strcmp(varargin{k}, 'x0')
                    self.x0 = varargin{k+1};
                elseif strcmp(varargin{k}, 'y0')
                    self.y0 = varargin{k+1};
                elseif strcmp(varargin{k}, 'shouldUseRobustFit')
                    self.shouldUseRobustFit = varargin{k+1};
                end
            end
            
            [self.x,self.y] = LatitudeLongitudeToTransverseMercator( self.lat, self.lon, self.lon0, self.k0 );
            if isnan(self.x0)
                self.x0 = min(self.x);
            end
            if isnan(self.y0)
                self.y0 = min(self.y);
            end
            
            self.x = self.x - self.x0;
            self.y = self.y - self.y0;
            
            self.xbar = mean(self.x);
            self.ybar = mean(self.y);
            xtilde = (self.x - self.xbar);
            ytilde = (self.y - self.ybar);
            M_xx = mean(xtilde.*xtilde);
            M_yy = mean(ytilde.*ytilde);
            M_xy = mean(xtilde.*ytilde);
            [A, ~] = eig([M_xx, M_xy; M_xy, M_yy]);
            v = A(:,2);
            self.theta = atan(v(2)/v(1));% + pi/4;
            self.q = xtilde*cos(self.theta) + ytilde*sin(self.theta);
            self.r = -xtilde*sin(self.theta) + ytilde*cos(self.theta);
            
            self.distribution = StudentTDistribution(self.sigma,self.nu);
            % empirically determined autocorrelation function
            self.distribution.rho = @(dt) exp(max(-abs(dt)/100., - abs(dt)/760 -1.3415)) ;
            
            if self.shouldUseRobustFit == 1
                self.spline_x = RobustTensionSpline(self.t,self.q,self.distribution,'K',self.K,'T',self.T,'lambda',Lambda.fullTensionIterated);
                self.spline_y = RobustTensionSpline(self.t,self.r,self.distribution,'K',self.K,'T',self.T,'lambda',Lambda.fullTensionIterated);
            else
                self.spline_x = TensionSpline(self.t,self.q,self.distribution,'K',self.K,'T',self.T);
                self.spline_y = TensionSpline(self.t,self.r,self.distribution,'K',self.K,'T',self.T);
            end
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % epsilon = observed position minus the fit position
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function epsilon = epsilon(self)
            epsilon = cat(2,self.spline_x.epsilon,self.spline_y.epsilon);
        end
        
        function a = isConstrained(self)
           a = self.spline_x.isConstrained && self.spline_y.isConstrained;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Projection
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [x,y] = xyAtTime(self,time)
            q_ = self.spline_x(time);
            r_ = self.spline_y(time);
            
            x = q_*cos(self.theta) - r_*sin(self.theta) + self.xbar;
            y = q_*sin(self.theta) + r_*cos(self.theta) + self.ybar;
        end
        
        function [lat,lon] = latitudeLongitudeAtTime(self,time)
            [x_,y_] = self.xyAtTime(time);
            [lat,lon] = TransverseMercatorToLatitudeLongitude(x_+self.x0,y_+self.y0,self.lon0, self.k0);
        end
        
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
        % 2D methods for outlier identification
        %
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [outlierDistribution,alpha] = setSigmaFromOutlierDistribution(self,rejectionPDFRatio)
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
            epsilon_full = cat(1,self.spline_x.epsilon,self.spline_y.epsilon);
            [outlierDistribution,alpha] = RobustTensionSpline.estimateOutlierDistributionFromKnownNoise(epsilon_full,self.distribution);
%             if self.alpha > 0
%                 self.sigma = RobustTensionSpline.sigmaFromOutlierDistribution(self.alpha,self.outlierDistribution,self.noiseDistribution,epsilon_full,rejectionPDFRatio);
%             end
        end
        
        function tensionValue = minimize(self,penaltyFunction)
            % the penalty function *must* take a tension spline object and
            % return a scalar. The function will be minimized by varying
            % lambda.
            [~,tensionValueBelow,tensionValueAbove] = GPSTensionSpline.minimizeFunctionOfSplineWithGridSearch(self,penaltyFunction);
            tensionValue = mean([tensionValueAbove,tensionValueBelow]);
            self.tensionValue = tensionValue;
%             tensionValue = GPSTensionSpline.minimizeFunctionOfSplineBounded(self,penaltyFunction,tensionValueBelow,tensionValueAbove);
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
                    
                    varargout{1} = self.spline_x(time,NumDerivatives);
                    varargout{2} = self.spline_y(time,NumDerivatives);
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GET %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                case '.'
                    [varargout{1:nargout}] = builtin('subsref',self,index);
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RESTRICT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                case '{}'
                    error('The GPSTensionSpline class does not know what to do with {}.');
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
        %
        % This algorithm also *stops* looking at higher lambdas once the
        % problem becomes constrained.
        function [lambda,lambdaBelow,lambdaAbove] = minimizeFunctionOfSplineWithGridSearch(aBivariateSpline,functionOfBivariateSpline,dlog10tensionValue,n)
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
            errorFunction = @(log10tensionValuePlusEpsilon) GPSTensionSpline.functionOfSplineWrapper(aBivariateSpline,log10tensionValuePlusEpsilon,functionOfBivariateSpline);
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
