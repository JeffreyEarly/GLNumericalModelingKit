function matrix = FourierTransformForwardMatrixReal(N)
% CosineTransformForwardMatrix  Discrete Cosine Transform (DCT-I) matrix
%
% This matrix exactly matches CosineTransformForward. See its documentation
% for details.

matrix = zeros(N/2+1,N);

for k=1:(N/2+1)
	for j=1:N
		matrix(k,j) = (1/N)*exp(-sqrt(-1)*2*pi*(j-1)*(k-1)/N);	
	end
end

return