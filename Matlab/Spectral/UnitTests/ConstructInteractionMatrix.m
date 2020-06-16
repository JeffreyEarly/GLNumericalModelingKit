function C = ConstructInteractionMatrix(F,G,H)
N = size(F,2);

C = zeros(N,N,N);
for i=1:N
    Fi = F(:,i);
    for j=1:N
        Gj = G(:,j);
        C(i,j,:) = H\(Fi.*Gj);
    end
end

end