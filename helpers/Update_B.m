function B = Update_B(FC, X, B_prev, M, R, xi_sq, Gam)

% permute FC to S*d*d
FC = permute(FC,[3 1 2]);


[~,q] = size(X);
[~,K,~] = size(B_prev);

B = B_prev;

ind_list = cell(K,1);

for k = 1:K
    ind_list{k} = find(M(k,:)~=0);
end

for k1 = 1:K
    
    for k2 = 1:K
        
        if k1 < k2 
            continue
        end
        
        i_list = ind_list{k1};
        j_list = ind_list{k2};
        
        if isempty(i_list) || isempty(j_list)
            continue
        end
        
        
        tmp_xi_sq = xi_sq(k1,k2,:);
        tmp_xi_sq = reshape(tmp_xi_sq,[q 1]);
        
        U = diag(1./tmp_xi_sq);
        V = zeros(1,q);
        
        tmp_empty = 0;
        for i = i_list
            for j = j_list
                
                ind1 = max(i, j);
                ind2 = min(i, j);
                
                if Gam(ind1,ind2) == 1                    
                    tmp_empty = 1;
                    
                    U = U + (X./R(:,ind1,ind2))'*X;
                    V = V + (FC(:,ind1,ind2)./R(:,ind1,ind2))'*X;
                end
                
            end
        end
        
        if tmp_empty == 1
            inv_U = pinv(U);
        
            tmp_b = mvnrnd(inv_U*(V.'), (inv_U + inv_U.')/2);

            B(k1,k2,:) = tmp_b;
        end
        
    end
end

% symmetric
for l = 1:q
    tmp_B = B(:,:,l);
    B(:,:,l) = tril(tmp_B) + tril(tmp_B,-1)';
end



end