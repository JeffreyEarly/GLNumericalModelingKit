function matrix = FourierTransformBackMatrix(N)
% CosineTransformBackMatrix  Discrete Cosine Transform (DCT-I) matrix
%
% This matrix exactly matches CosineTransformBack. See its documentation
% for details.

matrix = zeros(N,N);

for k=1:N
	for j=1:N
		matrix(k,j) = exp(sqrt(-1)*2*pi*(j-1)*(k-1)/N);	
	end
end

% matrix(:,1) = matrix(:,1)/2;

return