% Copyright (c) 2007, Stephen Morris 
% All rights reserved.
function roots = FindRootsFromChebyshevVector(f_cheb, zLobatto)
n = length(f_cheb);

A=zeros(n-1);   % "n-1" because Boyd's notation includes the zero-indexed
A(1,2)=1;       % elements whereas Matlab's of course does not allow this.
% In other words, to replicate Boyd's N=5 matrix we need to
% set n=6.
for j=2:n-2
    for k=1:n-1
        if j==k+1 || j==k-1
            A(j,k)=0.5;
        end
    end
end
for k=1:n-1
    A(n-1,k)=-f_cheb(k)/(2*f_cheb(n));  % c(1) in our notation is c(0) in Boyd's
end
A(n-1,n-2)=A(n-1,n-2)+0.5;
% Now we have the companion matrix, we can find its eigenvalues using the
% MATLAB built-in function.
eigvals=eig(A);

% We're only interested in the real elements of the matrix:
realvals=(arrayfun(@(x) ~any(imag(x)),eigvals)).*eigvals;

% Of course these are the roots scaled to the canonical interval [-1,1]. We
% need to map them back onto the interval [a,b]; we widen the interval just
% a fraction to make sure that we don't miss any that are right on the
% edge.
a = min(zLobatto);
b = max(zLobatto);
rangevals=nonzeros((arrayfun(@(x) abs(x)<=1.001, realvals)).*realvals);
roots=sort((rangevals.*0.5*(b-a)) + (0.5*(b+a)));
end