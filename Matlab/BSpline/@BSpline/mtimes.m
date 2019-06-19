function f = mtimes(f,g)
%.* BSpline multiplication

if ( ~isa(f, 'BSpline') )
    % Ensure BSpline is the first input:
    f = mtimes(g, f);
elseif ( isempty(g) )          % BSpline * []
    f = [];
elseif ( isnumeric(g) && isscalar(g) )
    h = BSpline(f.K,f.t_knot,f.m);
    h.x_std = g*f.x_std;
    h.x_mean = g*f.x_mean;
    f = h;
else
    error('This case is not handled!')
end
