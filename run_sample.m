
%% Load the sample data
sample_data_name = 'sample_data.mat';
data = load_data(sample_data_name, {'FC', 'age', 'sex', 'site_label'});

%% specify required parameters
option.thres = 10; % Threshold for edge screening (e.g., 10 or 5). 
                   % The screening rule uses 'R-squared > thres%'.
option.mite = 2000; % Total MCMC iterations
option.K = 'auto'; % Inital number of clusters, can be 'auto' or a numerical value. 
                   % If K is set to 'auto', then K is determined by the bulit-in method. 

option.output = 'outputs3'; % Path for outputs     
            
%% Run MCMC algorithm for clustering-enabled regression
rng(1)

MCMC_outputs = clustering_enabled_reg(data, option);

%% Make infenrence on MCMC samples
posterior_outputs = posterior_inf(MCMC_outputs, option.output);

%% Plot 95% CI for population-mean trajectory between two region clusters
disp('---------- Visualization example ----------')
disp('Load data information')
load(strcat(option.output,'/preprocessed_data.mat'), 'preprocessed_data')

data_info = preprocessed_data.data_info;


disp('Load estimated coefficients for trajectories between region clusters')
load(strcat(option.output,'/posterior_inference_outputs.mat'), 'posterior_outputs')


B_all_sorted = posterior_outputs.B_all_sorted;
clusters = posterior_outputs.clusters;


disp('Trajectory plots between each pair of regions')
for k1 = 1:length(clusters)
    for k2 = k1:length(clusters)
        plot_traj(k1, k2, B_all_sorted, data_info, option.output);
    end
end


disp('---------- Done ----------')


