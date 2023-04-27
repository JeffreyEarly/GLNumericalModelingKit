classdef InterpolatingSpline < BSpline
    %InterpolatingSpline Summary of this class goes here
    %   2 argument initialization
    %       f = InterpolatingSpline(t,x);
    %   where
    %       t       array of values for the independent axis
    %       x       array of values for the dependent axis 
    %       f       cubic spline interpolant
    % 
    %   3 argument initialization
    %       f = InterpolatingSpline(t,x,K);
    %   where
    %       K       order of the spline
    
    
    properties (Access = public)

    end
    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Initialization
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function self = InterpolatingSpline(t,x,varargin)
            % Make sure (t,x) are the right shapes, and see how many
            % dimensions of x we are given.
            t = reshape(t,[],1);
            N = length(t);
            x = reshape(x,[],1);
            if length(x) ~= N
                error('x and t must have the same length.');
            end
            
            if any(isnan(x)) || any(isnan(t))
                error('Some of the data contain NaNs!')
            end
            
            % Parse any extra input options.
            nargs = length(varargin);
            if nargs == 1
                K = varargin{1};
            elseif nargs == 0
                K = 4;
            elseif mod(nargs,2) == 0
                for k = 1:2:length(varargin)
                    if strcmp(varargin{k}, 'K')
                        K = varargin{k+1};
                    elseif strcmp(varargin{k}, 'S')
                        K = varargin{k+1}+1;
                    end
                end
            else
                error('Arguments must be given as name/value pairs');
            end
            
            % create knot points for a canonical interpolating spline
            t_knot = InterpolatingSpline.KnotPointsForPoints(t,K);
            
            x_mean = mean(x);
            x = x - x_mean;
            
            x_std = std(x);
            if x_std > 0
                x = x/x_std;
            else
                x_std = 1;
            end
            
            % Find the spline coefficients
            X = BSpline.Spline( t, t_knot, K, 0 );
            m = X\x;
            
            self@BSpline(K,t_knot,m);
            self.x_mean = x_mean;
            self.x_std = x_std;
        end

    end
    
    
    methods (Static)
                
        function t_knot = KnotPointsForDataPoints( t, options)
            arguments
                t (:,1) double
                options.K (1,1) double {mustBePositive,mustBeInteger,mustBeGreaterThanOrEqual(options.K,1)} = 4
                options.M (1,1) double {mustBePositive,mustBeInteger} = length(t)
            end
            mustBeGreaterThanOrEqual(options.M,options.K);
            mustBeLessThanOrEqual(options.M,length(t));

            N = length(t);
            t_pseudo = interp1((0:N-1)',t,linspace(0,N-1,options.M).');
            K = options.K;

            if mod(K,2) == 1
                % Odd spline order, so knots go in between points.
                dt = diff(t_pseudo);

                % This gives us N+1 knot points
                t_knot = [t_pseudo(1); t_pseudo(1:end-1)+dt/2; t_pseudo(end)];

                % Now remove start and end knots
                for i=1:((K-1)/2)
                    t_knot(2) = [];
                    t_knot(end-1) = [];
                end

            else
                t_knot = t_pseudo;

                % Now remove start and end knots
                for i=1:((K-2)/2)
                    t_knot(2) = [];
                    t_knot(end-1) = [];
                end

            end
        
            % Now we increase the multiplicity of the knot points at the beginning and
            % the end of the interval so that the splines do not extend past the end
            % points.
            t_knot = [repmat(t_knot(1),K-1,1); t_knot; repmat(t_knot(end),K-1,1)];
        end

        function t_knot = KnotPointsForPoints( t, K, DF )
            %% KnotPointsForPoints
            %
            % Returns the natural knot points for splines of order K for
            % observations at time points t. These knot points are place so
            % that each spline has exactly enough support for a fully
            % constrained system.
            %
            % It defaults to 1 degree of freedom, but you can override this
            % by setting DF to some nonnegative integer value.
            %
            % This matches the definitions in interpolation with tension
            % splines paper.
            
            if nargin < 3
                DF = 1;
            end
            
            if (DF < 1 || mod(DF,1) ~= 0)
                disp('DF must be a non-negative integer');
                return;
            end
            
            t = sort(t);
            
            t = [t(1); t(1+DF:DF:end-DF); t(end)];
            
            if length(t) > 2
                if mod(K,2) == 1
                    % Odd spline order, so knots go in between points.
                    dt = diff(t);

                    % This gives us N+1 knot points
                    t_knot = [t(1); t(1:end-1)+dt/2; t(end)];

                    % Now remove start and end knots
                    for i=1:((K-1)/2)
                        t_knot(2) = [];
                        t_knot(end-1) = [];
                    end

                else
                    t_knot = t;

                    % Now remove start and end knots
                    for i=1:((K-2)/2)
                        t_knot(2) = [];
                        t_knot(end-1) = [];
                    end

                end
            else
                t_knot = t;
            end
            
            % Now we increase the multiplicity of the knot points at the beginning and
            % the end of the interval so that the splines do not extend past the end
            % points.
            t_knot = [repmat(t_knot(1),K-1,1); t_knot; repmat(t_knot(end),K-1,1)];
        end
        
        function t_knot = KnotPointsForSplines( t, K, N_splines )
            %% KnotPointsForSplines
            %
            % Creates knot points for the data such that there will be
            % N_splines, each of which have approximately the same number
            % of observation points in support.
            %
            % Very crude, quick way of doing this.
            
            if N_splines < K
                N_splines = K;
            end
            
            t_knot = InterpolatingSpline.KnotPointsForPoints(t, K, ceil(length(t)/N_splines));
        end
  
    end
end


