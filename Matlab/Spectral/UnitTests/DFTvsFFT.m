% fprintf('Testing Fourier transform:\n')
% N=32; % total points
% T=4.0; % total time length
% t=T*(0:(N-1))'/N;
% 
% % Create a simple function, x = 2 + 3*sin(4*pi*t);
% % The lowest resolved frequency is one cycle per unit time
% x = 1*ones(size(t)) + 3*sin(2*(2*pi)*t);
% 
% [xbar, f] = FourierTransformForward(t,x);
% 
% S = T*(xbar .* conj(xbar));
% 
% dt = t(2)-t(1);
% df = f(2)-f(1);
% 
% x_sum = (1/T)*sum(x.*x)*dt;
% S_sum = sum(S)*df;
% 
% if abs(x_sum - S_sum) < 1e-7
%    fprintf('\tPower matches.\n')
% else
%    fprintf('\tPower FAIL\n')
% end
% 
% [x_back, t_back] = FourierTransformBack(f,xbar,1);
% if (max(abs(x-x_back)) < 1e-7) && (max(abs(t-t_back)) < 1e-7)
%    fprintf('\tTransformed backed.\n')
% else
%    fprintf('\tInverse transform FAIL.\n')
% end

N = 256;

DFT = fft(eye(N)).*1/N;
iDFT = N*ifft(eye(N));

DFTrp = real(DFT);
DFTip = imag(DFT);

rDFT = DFT(1:(N/2 + 1 ),:);
rDFTrp = real(rDFT);
rDFTip = imag(rDFT);

irDFT = iDFT(:,1:(N/2 + 1 ));
irDFT(:,2:(N/2)) = 2*irDFT(:,2:(N/2));

irDFTrp = real(irDFT);
irDFTip = -imag(irDFT);



NyNz = 128*128;
x = rand(N,NyNz);
nLoops = 10;

tic
for i=1:nLoops
    xbar_rp = rDFTrp*x;
    xbar_ip = rDFTip*x;
    xback = irDFTrp*xbar_rp + irDFTip*xbar_ip;
end
toc

tic
for i=1:nLoops
    xbar0 = rDFT*x;
    xback0 = irDFTrp*real(xbar0) + irDFTip*imag(xbar0);
end
toc

tic
for i=1:nLoops
    xbar1 = rDFT*x;
    xback1 = real(irDFT*xbar1);
end
toc

% tic
% xbar_rp = DFTrp*x;
% xbar_ip = DFTip*x;
% xbar0 = xbar_rp + sqrt(-1)*xbar_ip;
% toc

tic
for i=1:nLoops
    xbar2 = DFT*x;
    xback2 = iDFT*xbar2;
end
toc

tic
for i=1:nLoops
    xbar3 = fft(x)/N;
    xback3 = ifft(xbar3)*N;
end
toc

% xback = irDFTrp*real(xbar) + irDFTip*imag(xbar);