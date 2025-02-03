function [B_all, M_all, p_all, sigma_sq_all, delta_sq_all] = ...
    Gibbs_sampling(mite, FC, X, B0, M0, p0, sigma_sq0, delta_sq0,...
    Gam0, ind_batch, alpha, rho0, varrho0)

% Store Initialized Parameters
[d,~,S] = size(FC);
[~,q] = size(X);
[K,~] = size(M0);
[~,num_batch] = size(ind_batch);

ite = 1;

B_all = NaN(K,K,q,mite+1); 
B_all(:,:,:,ite) = B0;

M_all = NaN(K,d,mite+1);
M_all(:,:,ite) = M0;

p_all = NaN(K,mite+1);
p_all(:,ite) = p0;

sigma_sq_all = NaN(d,d,mite+1);
sigma_sq_all(:,:,ite) = sigma_sq0;

delta_sq_all = NaN(d,d,num_batch-1,mite+1);
delta_sq_all(:,:,:,ite) = delta_sq0;

tmp_FC = permute(FC,[3 1 2]); 

%% Run Gibbs Sampling
while ite <= mite
    if mod(ite, 10) == 0
        fprintf('MCMC Iteration %d \n', ite)
    end
    
    %% Update R
    %disp('Update R')
    [tmp_R0, tmp_R] = ...
        Update_R(sigma_sq_all(:,:,ite), delta_sq_all(:,:,:,ite),...
        Gam0, ind_batch);
    
    
    
    %% Update M
    %disp('Update M')
    tmp_M_all = M_all(:,:,ite);

    i_permu_list = randperm(d);

    for ii = 1:d
        
        i = i_permu_list(ii);
        
        tmp_FC_i = NaN(S,d);
        tmp_Gam_i = NaN(d,1);
        for j = 1:d
            if i == j 
                continue
            end
            ind1 = max(i, j);
            ind2 = min(i, j);
            
            tmp_FC_i(:,j) = tmp_FC(:,ind1,ind2);
            tmp_Gam_i(j) = Gam0(ind1,ind2);
        end
                
        tmp_R_i = tmp_R(:,:,i); 

        tmp_M_all(:,i) = ...
            Update_m(tmp_FC_i, X, B_all(:,:,:,ite), tmp_M_all, p_all(:,ite),...
            tmp_R_i, tmp_Gam_i);
        
    end
    M_all(:,:,ite+1) = tmp_M_all;
    
    %% Update xi_sq
    %disp('Update xi_sq')
    tmp_xi_sq_all = Update_xi_sq(B_all(:,:,:,ite), rho0);
    
    
    %% Update B
    %disp('Update B')
    B_all(:,:,:,ite+1) =...
        Update_B(FC, X, B_all(:,:,:,ite), M_all(:,:,ite+1),...
        tmp_R, tmp_xi_sq_all, Gam0);
        
    %% Update p
    %disp('Update p')
    p_all(:,ite+1) = Update_p(M_all(:,:,ite+1), alpha);

    
    %% Update sigma_sq
    %disp('Update sigma_sq')
    sigma_sq_all(:,:,ite+1) = ...
        Update_sigma_sq(FC, X, B_all(:,:,:,ite+1), M_all(:,:,ite+1), tmp_R0,...
        Gam0, varrho0);
    
    %% Update delta_sq
    %disp('Update delta_sq')
    delta_sq_all(:,:,:,ite+1) = ...
        Update_delta_sq(FC, X, B_all(:,:,:,ite+1), M_all(:,:,ite+1), sigma_sq_all(:,:,ite+1),...
        Gam0, ind_batch, varrho0);

    %%
    ite = ite + 1;
    
end




end