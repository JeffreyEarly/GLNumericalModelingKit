function [a] = ChebyshevCoefficientsFromIntegration(x, f, N )
%% Chebyshev Coefficients from Integration
% Compute the first N Chebyshev coefficients of the function f defined on
% grid x by integrating.
    x_norm = ChebyshevPolynomialsOnGrid( x );
    
    a = zeros(N,1);
    for i=1:N
        % This takes the interpolated function and multiplies it by the
        % (i-1)th Chebyshev Polynomial.
        %fun = @(y) interp1(x_norm,f,y,'spline').*cos((i-1)*acos(y))./sqrt(1-y.*y);
        
        % This takes the interpolated function, subtracts off the truncated
        % Chebyshev series, then multiplies by the next Chebyshev
        % Polynomial. In contrast to the previous method, this has much
        % better convergence properties for the integral.
        fun = @(y) (interp1(x_norm,f,y,'spline')-sum(repmat(a(1:i),1,length(y)).*cos(((0:i-1)')*acos(y)),1)).*cos((i-1)*acos(y))./sqrt(1-y.*y);
        %a(i) = quadgk(fun,-1,1,'MaxIntervalCount',10000);
        a(i) = integral(fun,-1,1)/(pi/2);
        if i==1
            a(1)=a(1)/2;
        end
    end
%     a = a/(pi/2);
%     a(1) = a(1)/2;
end
