function [zLobatto, f_cheb] = ProjectOntoChebyshevPolynomialsWithTolerance(f, bounds, tol)
m = 3;
m_max = 15;

Lz = max(bounds)-min(bounds);
zMin = min(bounds);

n = 2^m + 1;
cutoff = n;
while (cutoff == n && m <= m_max)
    m = m + 1;
    n = 2^m + 1;
    
    zLobatto = (Lz/2)*( cos(((0:n-1)')*pi/(n-1)) + 1) + zMin;
    f_cheb = fct(f(zLobatto));
    cutoff = standardChop(f_cheb, tol);
end

if cutoff < n
    n = cutoff;
    zLobatto = (Lz/2)*( cos(((0:n-1)')*pi/(n-1)) + 1) + zMin;
else
    disp('Unable to project density function within requested tolerance! Using maximum allowed length.');
end

f_cheb = f_cheb(1:n);
end