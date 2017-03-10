function u = ifct(uh)
%% Inverse Fast Chebyshev Transform 

N = length(uh) - 1;
s = N*[uh(1)*2; uh(2:N); uh(end)*2];
u = ifft([s; flip(s(2:end-1))],'symmetric');

u=u(1:N+1);

return
