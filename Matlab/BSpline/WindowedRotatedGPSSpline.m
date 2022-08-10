classdef WindowedRotatedGPSSpline < handle
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
        t
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
        splines_x = {}
        splines_y = {}
        pointSplineIndex
        splineParams = {}
    end
        
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Initialization
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function self = WindowedRotatedGPSSpline(t_raw,lat,lon,options)
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
            self.t = t;
               
            [x,y] = LatitudeLongitudeToTransverseMercator( self.lat, self.lon, self.lon0, self.k0 );  
            x = x - self.x0;
            y = y - self.y0;


            % ideally length(t)/(windowLength/3)
            nWindows = max(floor(length(t)/self.minOverlap - 2),1);
            startIndex = 1;
            self.pointSplineIndex = zeros(size(t));
            for iWindow = 1:nWindows
                endIndex = startIndex + self.minOverlap*3;
                if iWindow == nWindows
                    endIndex = length(self.t);
                end
                indices = startIndex:endIndex;
                [splineX,splineY,theta,x1,y1] = self.rotatedSplines(self.t(indices),x(indices),y(indices));
                self.splines_x{iWindow} = splineX;
                self.splines_y{iWindow} = splineY;
                self.splineParams{iWindow} = struct('theta',theta,'x1',x1,'y1',y1);
                if nWindows == 1
                    self.pointSplineIndex(:) = 1;
                elseif nWindows == 2
                    if iWindow == 1
                        self.pointSplineIndex(1:floor(length(t)/2)) = 1;
                    else
                        self.pointSplineIndex((floor(length(t)/2)+1):end) = 2;
                    end
                else
                    if iWindow == 1
                        self.pointSplineIndex(startIndex:(startIndex+2*self.minOverlap-1)) = iWindow;
                    elseif iWindow == nWindows
                        self.pointSplineIndex((startIndex+self.minOverlap):end) = iWindow;
                    else
                        self.pointSplineIndex((startIndex+self.minOverlap):(startIndex+2*self.minOverlap-1)) = iWindow;
                    end
                end
                startIndex = startIndex + self.minOverlap;
            end

        end

        function [splineX,splineY,theta,x1,y1] = rotatedSplines(self,t,x,y)
            % Fit a constant velocity through the points
            K = 2;
            t_knot = cat(1,min(t)*ones(K,1),max(t)*ones(K,1));
            splineMeanX = ConstrainedSpline(t,x,K,t_knot,self.gpsNoiseDistribution,[]);
            splineMeanY = ConstrainedSpline(t,y,K,t_knot,self.gpsNoiseDistribution,[]);
            dx = diff(splineMeanX([min(t) max(t)]));
            dy = diff(splineMeanY([min(t) max(t)]));
            x1 = splineMeanX(min(t));
            y1 = splineMeanY(min(t));
            theta = atan2(dy,dx);

            % rotate our points to align the mean velocity in x
            xp = (x-x1)*cos(theta) + (y-y1)*sin(theta);
            yp = -(x-x1)*sin(theta) + (y-y1)*cos(theta);
            u = sqrt(dx^2 + dy^2)/(max(t)-min(t));

            % splineX = SmoothingSpline(t,xp,self.gpsNoiseDistribution,'lambda',Lambda.fullTensionExpected);
            % splineY = SmoothingSpline(t,yp,self.gpsNoiseDistribution,'lambda',Lambda.fullTensionExpected);
            % alpha = 1/1000;
            % splineX.minimizeExpectedMeanSquareErrorInPercentileRange(alpha/2,1-alpha/2);
            % splineY.minimizeExpectedMeanSquareErrorInPercentileRange(alpha/2,1-alpha/2);


            splineY = SmoothingSpline(t,yp,self.gpsNoiseDistribution,'lambda',Lambda.fullTensionExpected);
            alpha = 1/1000;
            splineY.minimizeExpectedMeanSquareErrorInPercentileRange(alpha/2,1-alpha/2);

            splineX = SmoothingSpline(t,xp,self.gpsNoiseDistribution);
%             alpha = 1/10000000;
%             splineX.minimizeExpectedMeanSquareErrorInPercentileRange(alpha/2,1-alpha/2);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Evaluation
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [x,y] = xyAtTime(self,time,options)
            arguments
                self
                time (:,1) double
                options.NumDerivatives (1,1) double = 0
            end
            splineIndices = interp1(self.t,self.pointSplineIndex,time,'nearest');
            splineIndices(time>max(self.t)) = self.pointSplineIndex(end);
            splineIndices(time<min(self.t)) = self.pointSplineIndex(1);
            x = zeros(size(time));
            y = zeros(size(time));
            uniqueIndices = unique(splineIndices);
            for iUnique=1:length(uniqueIndices)
                splineIndex = uniqueIndices(iUnique);
                xp = self.splines_x{splineIndex}.ValueAtPoints(time(splineIndices==splineIndex), options.NumDerivatives);
                yp = self.splines_y{splineIndex}.ValueAtPoints(time(splineIndices==splineIndex), options.NumDerivatives);
                theta = -self.splineParams{splineIndex}.theta;
                xr = xp*cos(theta) + yp*sin(theta) + self.splineParams{splineIndex}.x1;
                yr = -xp*sin(theta) + yp*cos(theta) + self.splineParams{splineIndex}.y1;
                x(splineIndices==splineIndex) = xr;
                y(splineIndices==splineIndex) = yr;
            end
        end

        function [lat,lon] = latLonAtTime(self,time)
            [x,y] = self.xyAtTime(time);
            x = x + self.x0;
            y = y + self.y0;
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
                        [x,y] = self.xyAtTime(time,NumDerivatives=NumDerivatives);
                        varargout{1} = x;
                        varargout{2} = y;
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
