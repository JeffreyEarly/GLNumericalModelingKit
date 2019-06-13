function f = mtimes(f,g)
%.* BSpline multiplication

if ( ~isa(f, 'BSpline') )
    % Ensure BSpline is the first input:
    f = mtimes(g, f);
elseif ( isempty(g) )          % BSpline * []
    f = [];
elseif ( isnumeric(g) && isscalar(g) )
    f = BSpline(f.K,f.t_knot,g*f.m);
else
    error('This case is not handled!')
end
