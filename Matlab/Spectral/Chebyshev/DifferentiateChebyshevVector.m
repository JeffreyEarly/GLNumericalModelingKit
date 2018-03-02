function v_p = DifferentiateChebyshevVector(v)
v_p = zeros(size(v));
k = length(v)-1;
v_p(k) = 2*k*v(k+1);
for k=(length(v)-2):-1:1
    v_p(k) = 2*k*v(k+1) + v_p(k+2);
end
v_p(1) = v_p(1)/2;
end