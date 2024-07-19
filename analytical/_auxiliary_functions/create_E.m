theta = 0.25;

E = zeros(T,T);

for i = 1:T
    for j = 1:T
        if j < i
            E(i,j) = 1;
        else
            E(i,j) = theta^(j-i);
        end
    end
end