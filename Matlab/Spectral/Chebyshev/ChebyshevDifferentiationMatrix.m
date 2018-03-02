function D = ChebyshevDifferentiationMatrix(n)
%% Chebyshev Differentiation Matrix
% Returns the Chebyshev differentiation matrix for the first n polynomials.
    D = zeros(n,n);
    for i=1:n
        for j=1:n
            if ( j >= i+1 && mod(i+j,2)==1 )
                D(i,j) = 2*(j-1);
            else
                D(i,j) = 0.0;
            end
        end
    end
    D(1,:)=0.5*D(1,:);
end