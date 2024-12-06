N=129; % total points
T=8.0; % total time length
t=T*(0:(N-1))'/(N-1);

x = exp(-(t-T/2).^2);
figure, plot(t,x)

[x_s,f_s] = SineTransformForward(t,x);
[x_c,f_c] = CosineTransformForward(t,x);

t=T*(0:(N-1))'/N;
x = exp(-(t-T/2).^2);
[x_f,f_f] = FourierTransformForward(t,x);