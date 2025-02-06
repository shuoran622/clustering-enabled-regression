function initial_values = ...
    parameter_init(FC_norm, X_norm, ind_batch, Gam0, K, FC_test_pred_mat, output)

%%
disp('---------- Parameter initialization for Gibbs sampling ----------')

[d,~,~] = size(FC_norm);

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


%%
[d, ~] = size(Gam0);

Gam0_sym = tril(Gam0,-1) + tril(Gam0,-1)';

Gam0_sym_sum = zeros(d,1);
for i = 1:d
    Gam0_sym_sum(i) = nansum(Gam0_sym(i,:)) + nansum(Gam0_sym(:,i));
end

id_include_region = find(Gam0_sym_sum > 0);

Gam0 = Gam0(id_include_region,id_include_region);
sigma_sq0 = sigma_sq0(id_include_region,id_include_region);
delta_sq0 = delta_sq0(id_include_region,id_include_region,:);

FC_norm = FC_norm(id_include_region,id_include_region,:);

M0 = M0(:,id_include_region);

id_include_cluster = find(sum(M0,2) > 0);

K2 = length(id_include_cluster);

disp('Exclude regions with no edges')
disp(strcat('Remaining number of regions:', {' '}, string(size(M0,2))));

% If remove an entire cluster
if K2 < K  
    disp('')
    M0 = M0(id_include_cluster,:);
    B0 = B0(id_include_cluster,id_include_cluster,:);
    alpha = 1/K2.*ones(K2,1);
    p0 = alpha;
    disp(strcat('Number of clusters decreased due to excluding regions'))
    disp(strcat('Number of clusters:', {' '}, string(size(M0,1))));
else
    disp(strcat('Number of clusters:', {' '}, string(K)));
end


%%
initial_values = struct(...
    'FC_norm', FC_norm,...
    'B0', B0,...
    'sigma_sq0', sigma_sq0,...
    'M0', M0,...
    'Gam0', Gam0,...
    'delta_sq0', delta_sq0,...
    'p0', p0,...
    'alpha', alpha,...
    'id_include_region', id_include_region,...
    'id_include_cluster', id_include_cluster, ...
    'K', K, ...
    'd', d);
    


save_path = strcat(output,'/initial_values.mat');
save(save_path, 'initial_values', '-v7.3');


fprintf('Initial values saved as %s\n', save_path);


disp('---------- Done ----------')


end