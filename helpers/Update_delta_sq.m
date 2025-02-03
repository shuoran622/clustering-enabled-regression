function delta_sq = Update_delta_sq(FC,X,B,M,sigma_sq,Gam,ind_batch,varrho0)


% permute FC to S*d*d
FC = permute(FC,[3 1 2]);

[~,d,~] = size(FC);
[~,q] = size(X);

% tmp_b
tmp_b_summary = NaN(q,d,d);
for l = 1:q
    tmp_b_summary(l,:,:) = M'*B(:,:,l)*M; %K*d
end

[~,num_batch] = size(ind_batch);
delta_sq = NaN(d,d,num_batch-1);

for num_b = 1:(num_batch-1)
    
    tmp_ind_batch = ind_batch(:,num_b);

    %% batch indicators
    tmp_S = sum(tmp_ind_batch == 1);
    
    a = tmp_S/2 + varrho0;

    for i = 2:d

        for j = 1:(i-1)

            if Gam(i,j) == 1

                tmp = FC(:,i,j) - X*tmp_b_summary(:,i,j);

                tmp_batch = tmp(tmp_ind_batch == 1);

                b = 1/2.*tmp_batch.'*tmp_batch./sigma_sq(i,j) + varrho0;

                delta_sq(i,j,num_b) = 1 / gamrnd(a, 1/b);

            end

        end

    end

end

end