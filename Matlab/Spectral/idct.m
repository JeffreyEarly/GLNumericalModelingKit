function matrix = idct(n)

matrix = zeros(n,n);

for j=1:n
	for k=1:n
		matrix(k,j) = cos(pi*(j-1)*(k-1)/(n-1));	
	end
end

return