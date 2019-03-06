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
        
        shouldUseRobustFit = 0;
        
        distribution
        
        % splines
        N_max = 750 % maximum time series length
        N_fits
        t_boundary % boundaries defining which spline to use
        spline_x % cell array of TensionSplines, size(N_fits,1)
        spline_y % cell array of TensionSplines, size(N_fits,1)
        
        % t-distribution parameters
        % variance_of_the_noise = sigma*sigma*nu/(nu-2)
        sigma = 8.5
        nu = 4.5
        
        k0 = 0.9996
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
            
            self.distribution = StudentTDistribution(self.sigma,self.nu);
            % empirically determined autocorrelation function
%             self.distribution.rho = @(dt) exp(max(-abs(dt)/100., - abs(dt)/760 -1.3415)) ;
            
            minOverlapTime = 1200;
            minOverlapPoints = 2*self.K;
            minOverlap = max(round(minOverlapTime/median(diff(self.t))),minOverlapPoints);
            
            if 2*minOverlap>self.N_max
                self.N_max = 2*minOverlap;
                fprintf('Data is dense.\n',self.N_fits); 
            end
            
            self.N_fits = floor((N+minOverlap)/(self.N_max+minOverlap))+1;
            N_overlaps = self.N_fits-1;
            N_fit_length = floor((N - N_overlaps*minOverlap)/self.N_fits)+minOverlap;
            
            fit_starts = [1; 1+(1:(self.N_fits-1))'*N_fit_length-minOverlap] ;
            fit_ends = fit_starts+N_fit_length;
            fit_ends(end) = N;
                        
            if self.N_fits > 1
               fprintf('Splitting into %d overlapping fits.\n',self.N_fits); 
            end
            
            self.t_boundary = zeros(length(fit_starts)+1,1);
            self.t_boundary(1) = self.t(1);
            self.t_boundary(2:end-1) = self.t(fit_starts(2:end)) + (self.t(fit_ends(1:end-1))-self.t(fit_starts(2:end)))/2;
            self.t_boundary(end) = self.t(end);
            
            self.spline_x = cell(self.N_fits,1);
            self.spline_y = cell(self.N_fits,1);
            for iFit=1:self.N_fits
                indices = (fit_starts(iFit):fit_ends(iFit))';
                if self.shouldUseRobustFit == 1
                    self.spline_x{iFit} = RobustTensionSpline(self.t(indices),self.x(indices),self.distribution,'K',self.K,'T',self.T);
                    self.spline_y{iFit} = RobustTensionSpline(self.t(indices),self.y(indices),self.distribution,'K',self.K,'T',self.T);
                else
                    self.spline_x{iFit} = TensionSpline(self.t(indices),self.x(indices),self.distribution,'K',self.K,'T',self.T);
                    self.spline_y{iFit} = TensionSpline(self.t(indices),self.y(indices),self.distribution,'K',self.K,'T',self.T);
                end
            end
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Projection
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [x,y] = xyAtTime(self,time)
            x=zeros(size(time));
            y=zeros(size(time));
            spline_id = discretize(time,self.t_boundary);
            unique_spline_id = unique(spline_id);
            for iSpline=1:length(unique_spline_id)
                indices = (spline_id == unique_spline_id(iSpline));
                x(indices) = self.spline_x{unique_spline_id(iSpline)}(time(indices));
                y(indices) = self.spline_y{unique_spline_id(iSpline)}(time(indices));
            end
        end
        
        function [lat,lon] = latitudeLongitudeAtTime(self,time)
            [x_,y_] = self.xyAtTime(time);
            [lat,lon] = TransverseMercatorToLatitudeLongitude(x_+self.x0,y_+self.y0,self.lon0, self.k0);
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
            self.spline_x{1}.minimize(@(spline) self.distribution.andersonDarlingInterquartileError(spline.epsilon));
            self.spline_y{1}.minimize(@(spline) self.distribution.andersonDarlingInterquartileError(spline.epsilon));
            lastAlpha = 0.0;
            totalIterations = 0;
            [newOutlierDistribution, newAlpha] = RobustTensionSpline.estimateOutlierDistributionFromKnownNoise(cat(1,self.spline_x{1}.epsilon,self.spline_y{1}.epsilon),self.distribution);
            while (abs(lastAlpha-newAlpha) > 0.01 && totalIterations < 10)
                if newAlpha > 0 && ~isempty(newOutlierDistribution)
                    addedDistribution = AddedDistribution(newAlpha,newOutlierDistribution,self.distribution);
                else
                    addedDistribution = self.distribution;
                end
                self.spline_x{1}.minimize(@(spline) addedDistribution.andersonDarlingInterquartileError(spline.epsilon));
                self.spline_y{1}.minimize(@(spline) addedDistribution.andersonDarlingInterquartileError(spline.epsilon));
                lastAlpha = newAlpha;
                [newOutlierDistribution, newAlpha] = RobustTensionSpline.estimateOutlierDistributionFromKnownNoise(cat(1,self.spline_x{1}.epsilon,self.spline_y{1}.epsilon),self.distribution);
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
            epsilon_full = cat(1,self.spline_x{1}.epsilon,self.spline_y{1}.epsilon);
            [outlierDistribution,alpha] = RobustTensionSpline.estimateOutlierDistributionFromKnownNoise(epsilon_full,self.distribution);
%             if self.alpha > 0
%                 self.sigma = RobustTensionSpline.sigmaFromOutlierDistribution(self.alpha,self.outlierDistribution,self.noiseDistribution,epsilon_full,rejectionPDFRatio);
%             end
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
    
end
