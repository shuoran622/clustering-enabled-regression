function M_i = Update_m(FC_i, X, B, M_prev, p, R_i, Gam_i)


[~,d] = size(FC_i);
[~,q] = size(X);
[~,K,~] = size(B);

M = M_prev;


tmp_b_summary = NaN(q,K,d);
for l = 1:q
    tmp_b_summary(l,:,:) = B(:,:,l)*M; 
end


f_list = zeros(K,1);
    
for k = 1:K

    f_temp = 0;
    for j = 1:d
        if Gam_i(j) == 1

            FC_X_tmp_b = FC_i(:,j) - X*tmp_b_summary(:,k,j);
            
            f_temp = f_temp -1/2.*(FC_X_tmp_b./R_i(:,j))'*FC_X_tmp_b;

        end

    end

    f_list(k) = f_temp + log(p(k));

end

f_list = f_list - max(f_list).*ones(K,1);
f_list = exp(f_list);

f_list = f_list./sum(f_list);


M_i = mnrnd(1, f_list.', 1); 



end