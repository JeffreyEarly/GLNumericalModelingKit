% Want to create a chebyshev grid that never has two or more point between
% its points. If that makes sense.
function [z_lobatto_grid] = FindSmallestChebyshevGridWithNoGaps(z)
if (z(2) - z(1)) > 0 % make z_out decreasing
    z = flip(z);
end

L = max(z)-min(z);
np = ceil(log2(length(z)));
n = 2^np;
z_lobatto_grid = (L/2)*( cos(((0:n-1)')*pi/(n-1)) + 1) + min(z);

while( length(unique(interp1(z_lobatto_grid,z_lobatto_grid,z,'previous'))) ~= length(z) )
    np = np + 1;
    n = 2^np;
    z_lobatto_grid = (L/2)*( cos(((0:n-1)')*pi/(n-1)) + 1) + min(z);
end

end