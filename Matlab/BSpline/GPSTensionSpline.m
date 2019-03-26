classdef GPSTensionSpline < BivariateTensionSpline
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
        % bivariate gps data
        t_raw % format that was given
        lat
        lon
        
        % projected data
        lon0 % reference longitude of transverse mercator projection being used
        x0  % false easting
        y0  % false northing
        
        gpsNoiseDistribution = StudentTDistribution(8.5,4.5)
        
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
            t_raw = reshape(t,[],1);
            lat = reshape(lat,[],1);
            lon = reshape(lon,[],1);
            
            if isdatetime(t_raw)
                t = (datenum(t_raw) - datenum(t_raw(1)))*24*60*60;
            else
                t = t_raw;
            end
            
            if length(lat) ~= N || length(lon) ~= N
                error('lat, lon, and t must have the same length.');
            end
            
            lon0 = min(lon)+(max(lon)-min(lon))/2;
            x0 = 0;
            y0 = 0;
            gpsNoiseDistribution = StudentTDistribution(8.5,4.5);
            % works fine, but really slows it down
%             gpsNoiseDistribution.rho = @(dt) exp(max(-abs(dt)/100., - abs(dt)/760 -1.3415)) ;
            
            nargs = length(varargin);
            if mod(nargs,2) ~= 0
                error('Arguments must be given as name/value pairs');
            end
                        
            for k = 1:2:length(varargin)
                if strcmp(varargin{k}, 'gpsNoiseDistribution')
                    gpsNoiseDistribution = varargin{k+1};
                elseif strcmp(varargin{k}, 'lon0')
                    lon0 = varargin{k+1};
                elseif strcmp(varargin{k}, 'x0')
                    x0 = varargin{k+1};
                elseif strcmp(varargin{k}, 'y0')
                    y0 = varargin{k+1};
                end
            end
            
            k0 = 0.9996;
            [x,y] = LatitudeLongitudeToTransverseMercator( lat, lon, lon0, k0 );
            
            x = x - x0;
            y = y - y0;
            
            self@BivariateTensionSpline(t,x,y,gpsNoiseDistribution,varargin{:});
            self.gpsNoiseDistribution = gpsNoiseDistribution;
            self.t_raw = t_raw;
            self.lat = lat;
            self.lon = lon;
            self.lon0 = lon0;
            self.x0 = x0;
            self.y0 = y0;
            
           
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Evaluation
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [lat,lon] = latLonAtTime(self,time)
            x = self.spline_x(time) + self.spline_mean_x(time) + self.x0;
            y = self.spline_y(time) + self.spline_mean_y(time) + self.y0;
            [lat,lon] = TransverseMercatorToLatitudeLongitude(x,y,self.lon0,self.k0);
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
                    
                    if NumDerivatives == 0
                        [lat_,lon_] = self.latLonAtTime(time);
                        varargout{1} = lat_;
                        varargout{2} = lon_;
                    else
                        varargout{1} = self.spline_x(time,NumDerivatives) + self.spline_mean_x(time,NumDerivatives);
                        varargout{2} = self.spline_y(time,NumDerivatives) + self.spline_mean_y(time,NumDerivatives);
                    end
                    
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
