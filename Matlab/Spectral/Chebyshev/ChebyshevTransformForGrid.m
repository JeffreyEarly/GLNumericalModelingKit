% Given some Lobatto grid and some desired output grid, return the
% transformation function T that goes from spectral to the output
% grid. This basically gives you spectral interpolation.
function [T, doesOutputGridSpanDomain] = ChebyshevTransformForGrid(lobatto_grid, output_grid)
if (min(output_grid) == min(lobatto_grid) && max(output_grid) == max(lobatto_grid))
    doesOutputGridSpanDomain = 1;
else
    doesOutputGridSpanDomain = 0;
end

% T_out transforms vector solutions of the eigenvalue problem
% into gridded solution on z_out
if doesOutputGridSpanDomain == 1 && IsChebyshevGrid(output_grid) == 1
    if length(output_grid) == length(lobatto_grid)
        T = @(f_cheb) ifct(f_cheb);
    elseif length(output_grid) > length(lobatto_grid)
        T = @(f_cheb) ifct(cat(1,f_cheb,zeros(length(output_grid)-length(lobatto_grid),1)));
    elseif length(output_grid) < length(lobatto_grid)
        T = @(f_cheb) ifct(f_cheb(1:length(output_grid)));
    end
else
    L = max(lobatto_grid)-min(lobatto_grid);
    x = (2/L)*(output_grid-min(lobatto_grid)) - 1;
    t = acos(x);
    TT = zeros(length(output_grid),length(lobatto_grid));
    for iPoly=0:(length(lobatto_grid)-1)
        TT(:,iPoly+1) = cos(iPoly*t);
    end
    T = @(f_cheb) TT*f_cheb;
end
end