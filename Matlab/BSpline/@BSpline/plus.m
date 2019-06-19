function f = plus(f,g)
%+ BSpline addition

if ( ~isa(f, 'BSpline') )
    % Ensure BSpline is the first input:
    f = plus(g, f);
elseif ( isempty(g) )          % BSpline * []
    f = [];
elseif ( isnumeric(g) && isscalar(g) )
    %     X = f.B(:,:,1);
    %     m1 = X\ones(size(X,1),1);
    %     f = BSpline(f.K,f.t_knot,f.m + g*m1);
    h = BSpline(f.K,f.t_knot,f.m);
    h.x_std = f.x_std;
    h.x_mean = f.x_mean+g;
    f = h;
else
    error('This case is not handled!')
end
