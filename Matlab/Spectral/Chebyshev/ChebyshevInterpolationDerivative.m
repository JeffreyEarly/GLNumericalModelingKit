function D = ChebyshevInterpolationDerivative(n)
%% Chebyshev Interpolation Derivative
% taken from Canuto, et al. 2.4.33
    D = zeros(n,n);
    N = n-1;
    c = @(j) (j == 0 || j == N)*2 + (j>0 && j<N)*1;
    for j=0:(n-1)
        for l=0:(n-1)
            if j ~= l
               D(j+1,l+1) = -(c(j)/c(l))*((-1)^(j+l))/( sin( (j+l)*pi/(2*N) ) * sin( (j-l)*pi/(2*N) ) )/2;
            elseif j == l && j == 0
                D(j+1,l+1) = (2*N*N+1)/6;
            elseif j == l && j == N
                D(j+1,l+1) = -(2*N*N+1)/6;
            else
                D(j+1,l+1) = -cos(pi*j/N)/(2*(sin(j*pi/N))^2);
            end
        end
    end
    D(1,:)=D(1,:);
end