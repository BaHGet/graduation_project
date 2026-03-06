function vecs = linspace_vec(v1,v2,N)
vecs = zeros(N,2);
for i=1:N
    t = (i-1)/(N-1);
    vecs(i,:) = (1-t)*v1 + t*v2;
end
end