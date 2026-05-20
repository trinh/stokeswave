function M = henriciMatrix(m,n, delta)

M = zeros(m,n);

for k1=1:m
    for k2 = 1:n
        if k2 >= k1
            M(k1,k2) = nchoosek(k2-1,k1-1)*delta^(k2-k1);
        end
    end
end