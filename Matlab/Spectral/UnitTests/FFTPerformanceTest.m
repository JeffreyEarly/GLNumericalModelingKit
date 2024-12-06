N = 256;
NyNz = 256*256;
x = double(rand(N,NyNz));
nLoops = 400;

tic
for i=1:nLoops
    xbar3 = fft(x)/N;
    xback3 = ifft(xbar3,'symmetric')*N;
end
toc

% xback = irDFTrp*real(xbar) + irDFTip*imag(xbar);