function matrix = FourierTransformBackMatrixReal(N)
% FourierTransformBackMatrixReal  Discrete Fourier Transform matrix
%
% This matrix assumes you only have N/2+1 components and that the others
% are hermitian conjugate. You cannot directly apply the matrix and get the
% correct answer...
%
% rDFT = FourierTransformForwardMatrixReal(N);
% irDFT = FourierTransformBackMatrixReal(N);
%
% Option 1:
% xbar = rDFT*x;
% x = real(irDFT*xbar);
%
% Option 2:
% xbar = rDFT*x;
% x = irDFTrp*real(xbar) + irDFTip*imag(xbar);
% where
% irDFTrp = real(irDFT); irDFTip = -imag(irDFT);

matrix = zeros(N,N/2+1);

for k=1:N
	for j=1:(N/2+1)
		matrix(k,j) = exp(sqrt(-1)*2*pi*(j-1)*(k-1)/N);	
	end
end

matrix(:,2:(N/2)) = 2*matrix(:,2:(N/2));

% matrix(:,1) = matrix(:,1)/2;

return