
%% Load sample data
load('sample_data.mat')

%% specify required parameters
thres = 10; % thres/100 is for edge screening, e.g., 10, 5 
mite = 2000; % total MCMC iterations
K = 'auto'; % inital number of clusters, can be 'auto' or a numerical value. 
            % if K is set to 'auto', then K is determined by spectral clustering.
%% Run MCMC algorithm for clustering-enabled regression
rng(622)
[M_all, B_all, sigma_sq_all, delta_sq_all, edge_inclusion, data_info] = ...
    clustering_enabled_reg(FC, age, sex, site_label, thres, mite, K);

%% Make infenrence on MCMC samples
[clusters, mean_B, B_all_sorted, mean_sigma_sq, mean_delta_sq] = ...
    posterior_inf(M_all, B_all, sigma_sq_all, delta_sq_all);

%% Plot 95% CI for population-mean trajectory between two region clusters
k1 = 1;
k2 = 2;

plot_traj(k1, k2, B_all_sorted, data_info)