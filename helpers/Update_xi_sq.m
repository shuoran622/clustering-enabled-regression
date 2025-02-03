function xi_sq = Update_xi_sq(B, rho0)

[K,~,q] = size(B);

xi_sq = NaN(K,K,q);

a = 1/2 + rho0;

for l = 1:q
    
    for k1 = 1:K
        for k2 = 1:k1
            
            b = 1/2.*B(k1,k2,l)^2 + rho0;

            r = 1;
            while r > 0
                tmp = 1 / gamrnd(a, 1/b);

                r = double(tmp == Inf);
            end

            xi_sq(k1,k2,l) = tmp;
            
        end
    end
    
end


for l = 1:q
    tmp_xi_sq = xi_sq(:,:,l);
    xi_sq(:,:,l) = tril(tmp_xi_sq) + tril(tmp_xi_sq,-1)';
end


end