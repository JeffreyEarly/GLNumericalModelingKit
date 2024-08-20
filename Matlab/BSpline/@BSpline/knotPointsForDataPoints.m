function t_knot = knotPointsForDataPoints( t, options)
% create knot points appropriate for observation times t
%
% - Topic: Methodology (Static methods)
% - Declaration: t_knot = knotPointsForDataPoints( t, options)
% - Parameter t: observation times (N)
% - Parameter K: (optional) spline order
% - Parameter M: (optional) number of splines (M<=N)
% - Returns t_knot: vector of knot point locations
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