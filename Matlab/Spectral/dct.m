function matrix = dct(n)

matrix = zeros(n,n);

for k=1:n
	for j=1:n
		matrix(k,j) = (2/(n-1))*cos(pi*(j-1)*(k-1)/(n-1));	
	end
end

matrix(1,:) = matrix(1,:)/2;
matrix(n,:) = matrix(n,:)/2;
matrix(:,1) = matrix(:,1)/2;
matrix(:,n) = matrix(:,n)/2;

return