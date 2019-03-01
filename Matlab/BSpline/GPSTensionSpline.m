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
        t
        lat
        lon
        
        % projected data
        lon0 % reference longitude of transverse mercator projection being used
        x0  % false easting
        y0  % false northing
        x
        y
        
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
        
        k0 = 0.9996
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
            self.lat = reshape(lat,[],1);
            self.lon = reshape(lon,[],1);
            
            if length(x) ~= N || length(y) ~= N
                error('x, y, and t must have the same length.');
            end
            
            self.K = 4; % default spline order (cubic spline)
            self.T = self.K-1; % default tension *degree* (order-1)
            self.lon0 = min(lon)+(max(lon)-min(lon))/2;
            
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
                end
            end
            
            [self.x,self.y] = LatitudeLongitudeToTransverseMercator( self.lat, self.lon, self.lon0, self.k0 );
            self.x0 = min(self.x);
            self.y0 = min(self.y);
            self.x = self.x - self.x0;
            self.y = self.y - self.y0;
            
            self.spline_x = RobustTensionSpline(self.t,self.x,StudentTDistribution(self.sigma,self.nu));
            self.spline_y = RobustTensionSpline(self.t,self.y,StudentTDistribution(self.sigma,self.nu)); 
        end
        
        function [x,y] = xyAtTime(self,time)
            x = self.spline_x(time);
            y = self.spline_y(time);
        end
        
        function [lat,lon] = latitudeLongitudeAtTime(self,time)
            [x_,y_] = self.xyAtTime(time);
            [lat,lon] = TransverseMercatorToLatitudeLongitude(x_+self.x0,y_+self.y0,self.lon0, self.k0);
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
                    error('The GPSTensionSpline class does not know what to do with {}.');
                otherwise
                    error('Unexpected syntax');
            end
            
        end
    end
    
end
