classdef WindowedGPSSpline < handle
    % WindowedGPSSpline Fit noisy GPS data with a tensioned interpolating
    % spline
    %   3 argument initialization
    %       f = GPSSmoothingSpline(t,lat,lon);
    %   where
    %       t       time (seconds or datetime)
    %       lat     latitude
    %       lon     longitude
    %       f       spline interpolant
    %
    %   WindowedGPSSpline takes a number of optional input argument pairs.
    
    
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
        k0

        minOverlap
        spline_x
        spline_y
    end
        
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Initialization
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function self = WindowedGPSSpline(t_raw,lat,lon,options)
            arguments
                t_raw (:,1) double {mustBeNumeric,mustBeReal}
                lat (:,1) double {mustBeNumeric,mustBeReal,mustBeEqualSize(t_raw,lat)}
                lon (:,1) double {mustBeNumeric,mustBeReal,mustBeEqualSize(t_raw,lon)}
                options.gpsNoiseDistribution Distribution {mustBeNonempty} = StudentTDistribution(8.5,4.5)
                options.minOverlap = 10
                options.x0 (1,1) double {mustBeNumeric,mustBeReal} = 0
                options.y0 (1,1) double {mustBeNumeric,mustBeReal} = 0
                options.lon0 (1,1) double = min(lon)+(max(lon)-min(lon))/2;
                options.k0 (1,1) double = 0.9996;
            end
            self.t_raw = t_raw;
            self.lat = lat;
            self.lon = lon;
            self.lon0 = options.lon0;
            self.x0 = options.x0;
            self.y0 = options.y0;
            self.gpsNoiseDistribution = options.gpsNoiseDistribution;
            self.k0 = options.k0;
            self.minOverlap = options.minOverlap;
            
            if isdatetime(t_raw)
                t = (datenum(t_raw) - datenum(t_raw(1)))*24*60*60;
            else
                t = t_raw;
            end        
               
            [x,y] = LatitudeLongitudeToTransverseMercator( self.lat, self.lon, self.lon0, self.k0 );  
            x = x - self.x0;
            y = y - self.y0;
            self.spline_x = WindowedSmoothingSpline(t,x,self.gpsNoiseDistribution, minOverlap=self.minOverlap);
            self.spline_y = WindowedSmoothingSpline(t,y,self.gpsNoiseDistribution, minOverlap=self.minOverlap);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Evaluation
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [x,y] = xyAtTime(self,time)
            x = self.spline_x(time);
            y = self.spline_y(time);
        end

        function [lat,lon] = latLonAtTime(self,time)
            x = self.spline_x(time) + self.x0;
            y = self.spline_y(time) + self.y0;
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
                        varargout{1} = self.spline_x(time,NumDerivatives);
                        varargout{2} = self.spline_y(time,NumDerivatives);
                    end
                    
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
    
    
    
    
end

% Custom validation function
function mustBeEqualSize(a,b)
% Test for equal size
if ~isequal(size(a),size(b))
    eid = 'Size:notEqual';
    msg = 'Size of first input must equal size of second input.';
    throwAsCaller(MException(eid,msg))
end
end
