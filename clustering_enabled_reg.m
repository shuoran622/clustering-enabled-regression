function [M_all, B_all, sigma_sq_all, delta_sq_all, Gam0, data_info] = ...
    clustering_enabled_reg(FC, age, sex, site_label, thres, mite, K)

%%
% FC: funcitonal connectivity, a dxdxS array, d is the number of ROIs, S is the number of subjects.
% age is a Sx1 vector.
% sex is a Sx1 vector. 1 indicates males and 0 indicates females.
% site_label is a Sx1 vector, containing site label of corresponding fMRI scan. Label has to be integer starting from 1.
% thres: thres/100 is for edge screening, e.g., 10, 5 
% mite: total MCMC iterations
% K: inital number of clusters, can be 'auto' or a numerical value. If K is set to 'auto', then K is determined by spectral clustering.



%%


addpath('helpers'); 

% Construct X matrix from age and sex, and normalize FC and X. 
% Transfer site label into a matrix indicating batch indices.
[FC_norm, X_norm, ind_batch, data_info] = process_data(FC, age, sex, site_label);

% Gam0 is the matrix indicating the prescreen results,
% Gam0(i,j) = 1 indicates edge (i,j) has FC with R_squared above thres/100,
% thus included in the analysis
disp('Indepdent regression for edge screening...')
[Gam0, indepReg_info] = prescreen(FC_norm, X_norm, ind_batch, data_info, thres);
disp('---------done----------')

% Parameter Initialization
% Specify hyperparameters
rho0 = 10^-2;
varrho0 = 10^-6;
disp('Parameter initialization for gibbs sampling...')
[B0,sigma_sq0,M0,Gam0,delta_sq0,p0,K,alpha] = ...
    parameter_init(FC_norm, X_norm, ind_batch, Gam0, K, indepReg_info.FC_test_pred_mat);
disp('---------done----------')

% Exclude Regions with No Edges
disp('Excluding regions with no edges after edge screening...')

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

fprintf('Remaining number of regions is %d \n', size(M0,2))
% If remove an entire cluster
if K2 < K  
    M0 = M0(id_include_cluster,:);
    B0 = B0(id_include_cluster,id_include_cluster,:);
    alpha = 1/K2.*ones(K2,1);
    p0 = alpha;   
    fprintf('Number of clusters is adjusted to %d \n', size(M0,1))
end

disp('---------done----------')

%% MCMC
disp('Performing Gibbs Sampling...')

[B_all, M_all, ~, sigma_sq_all, delta_sq_all] = ...
    Gibbs_sampling(mite, FC_norm, X_norm, B0, M0, p0, sigma_sq0, delta_sq0,...
    Gam0, ind_batch, alpha, rho0, varrho0);

disp('---------done----------')


%% Adjust results based on id_include_region
d = size(FC,1);
if length(id_include_region) < d
    
    disp('Adjust results for excluded regions...')
    
    K = size(M_all,1);
    num_batch = size(ind_batch,2);
    tmp_M_all = zeros(K,d,mite+1); 
    tmp_M_all(1:K,id_include_region,:) = M_all;
    
    tmp_sigma_sq_all = NaN(d,d,mite+1);
    tmp_delta_sq_all = NaN(d,d,num_batch,mite+1);
    
    tmp_sigma_sq_all(id_include_region,id_include_region,:) = sigma_sq_all;
    tmp_delta_sq_all(id_include_region,id_include_region,1:(num_batch-1),:) = delta_sq_all;
    % add delta_sq = 1 for the last batch
    tmp_gam0 = Gam0;
    tmp_gam0(tmp_gam0 == 0) = NaN;
    for i = 1:(mite+1)
        tmp_delta_sq_all(id_include_region,id_include_region,num_batch,i) = tmp_gam0;
    end
    tmp_Gam0 = zeros(d,d);
    tmp_Gam0(id_include_region,id_include_region) = Gam0;
     
    M_all = tmp_M_all;   
    sigma_sq_all = tmp_sigma_sq_all;
    delta_sq_all = tmp_delta_sq_all;  
    Gam0 = tmp_Gam0;
    
    disp('---------done----------')
end


end

