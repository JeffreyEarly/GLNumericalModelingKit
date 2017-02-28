function bool = IsChebyshevGrid(z_in)
% make sure the grid is monotonically decreasing
if (z_in(2) - z_in(1)) > 0
    z_in = flip(z_in);
end

z_norm = ChebyshevPolynomialsOnGrid( z_in );
N_points = length(z_in);
xi=(0:N_points-1)';
lobatto_grid = cos(xi*pi/(N_points-1));
z_diff = z_norm-lobatto_grid;
if max(abs(z_diff)) < 1e-6
    bool = 1;
else
    bool = 0;
end
end