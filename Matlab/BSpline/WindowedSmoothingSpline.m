classdef WindowedSmoothingSpline < handle
% Choose a window length (the number of points that will be handed to a
% spline), and the middle 1/3 of that window (other than the ends) will be
% used to determined the smoothed values.
    properties
        t
        x
        minOverlap
        distribution
        splines = {};
        pointSplineIndex
    end

    methods

        function self = WindowedSmoothingSpline(t,x,distribution,options)
            arguments
                t (:,1) double {mustBeNonempty}
                x (:,1) double {mustBeNonempty}
                distribution Distribution {mustBeNonempty}
                options.minOverlap = 10
            end
            self.t = t;
            self.x = x;
            self.distribution = distribution;
            self.minOverlap = options.minOverlap;

            % ideally length(t)/(windowLength/3)
            nWindows = max(floor(length(t)/self.minOverlap - 2),1);
            startIndex = 1;
            self.pointSplineIndex = zeros(size(t));
            for iWindow = 1:nWindows
                endIndex = startIndex + self.minOverlap*3;
                if iWindow == nWindows
                    endIndex = length(t);
                end
                indices = startIndex:endIndex;
                self.splines{iWindow} = SmoothingSpline(t(indices),x(indices),distribution,'shouldUseRobustFit',1);
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
                        t_requested = idx{1};
                    end

                    if length(idx) >= 2
                        NumDerivatives = idx{2};
                    else
                        NumDerivatives = 0;
                    end
                    splineIndices = interp1(self.t,self.pointSplineIndex,t_requested,'nearest');
                    splineIndices(t_requested>max(self.t)) = self.pointSplineIndex(end);
                    splineIndices(t_requested<min(self.t)) = self.pointSplineIndex(1);
                    x_out = zeros(size(t_requested));
                    uniqueIndices = unique(splineIndices);
                    for iUnique=1:length(uniqueIndices)
                        splineIndex = uniqueIndices(iUnique);
                        x_out(splineIndices==splineIndex) = self.splines{splineIndex}.ValueAtPoints(t_requested(splineIndices==splineIndex), NumDerivatives);
                    end

                    varargout{1} = x_out;

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GET %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                case '.'
                    [varargout{1:nargout}] = builtin('subsref',self,index);

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RESTRICT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                case '{}'
                    error('The BSpline class does not know what to do with {}.');
                otherwise
                    error('Unexpected syntax');
            end

        end

    end

end