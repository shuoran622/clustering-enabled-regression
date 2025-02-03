function [B0,sigma_sq0,M0,Gam0,delta_sq0,p0,K,alpha] = ...
    parameter_init(FC_norm, X_norm, ind_batch, Gam0, K, FC_test_pred_mat)

[d,~,S] = size(FC_norm);

FC_permu = permute(FC_norm,[3 1 2]); % easy-to-use

% Cluster Initialization
[M0, K] = cluster_init(FC_test_pred_mat, K);

alpha = 1/K.*ones(K,1);
% Initializtion of p and delta squared
p0 = alpha;
delta_sq0 = ones(d,d,size(ind_batch,2)-1);

% Initialize B and sigma_sq by lm
q = size(X_norm,2);
B0 = zeros(K,K,q);
sigma_sq0 = NaN(d,d); 
tmp_e_list = NaN(K,K);

for k1 = 1:K
    for k2 = 1:K
        if k1 < k2 
            continue
        end
        
        i_list = find(M0(k1,:)~=0);
        j_list = find(M0(k2,:)~=0);
        
        if isempty(i_list) || isempty(j_list)
            continue
        end
        
        tmp_FC = [];
        tmp_X = [];
        
        for i = i_list
            for j = j_list
                
                ind1 = max(i,j);
                ind2 = min(i,j);
                
                if Gam0(ind1,ind2) == 1
                    
                    tmp_FC_ij = FC_permu(:,ind1,ind2);
                    
                    tmp_FC = [tmp_FC;tmp_FC_ij];
                    tmp_X = [tmp_X;X_norm];
                end

            end
        end
        
        if ~isempty(tmp_FC)
            
            tmp_md = fitlm(tmp_X,tmp_FC,'Intercept',false);
            
            tmptmp = tmp_md.Coefficients;

            tmp_pvalue = table2array(tmptmp(:,4));
            
            tmp_b = table2array(tmptmp(:,1)); 
            tmp_b(tmp_pvalue >= 0.05) = 0; 
            tmp_e = tmp_md.MSE;
            
            B0(k1,k2,:) = tmp_b;
            
            tmp_e_list(k1,k2) = tmp_e;
        end
        
    end
end


for k1 = 1:K      
    for k2 = 1:K
        
        if k1 < k2
            continue
        end
        
        i_list = find(M0(k1,:)~=0);
        j_list = find(M0(k2,:)~=0);
        
        
        for i = i_list
            for j = j_list
                
                ind1 = max(i,j);
                ind2 = min(i,j);
                
                if Gam0(ind1,ind2) == 1
                    sigma_sq0(ind1,ind2) = tmp_e_list(k1,k2);
                end
                
            end
        end
        
    end
end

for l = 1:q
    tmp_B = B0(:,:,l);
    B0(:,:,l) = tril(tmp_B) + tril(tmp_B,-1)';
end

B0(isnan(B0)) = 0;


end