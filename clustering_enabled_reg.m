function MCMC_outputs = ...
    clustering_enabled_reg(data, option)

%%
% Fields in data:
% FC: funcitonal connectivity, a dxdxS array, d is the number of ROIs, S is the number of subjects.
% age is a Sx1 vector.
% sex is a Sx1 vector. 1 indicates males and 0 indicates females.
% site_label is a Sx1 vector, containing site label of corresponding fMRI scan. Label has to be integer starting from 1.

% Fields in option: 
% thres: thres/100 is for edge screening, e.g., 10, 5 
% mite: total MCMC iterations
% K: inital number of clusters, can be 'auto' or a numerical value. If K is set to 'auto', then K is determined by spectral clustering.
% output: path for Path for outputs


%%
FC = data.FC;
age = data.age;
sex = data.sex;
site_label = data.site_label;

if ~isequal(unique(site_label), (1:length(unique(site_label)))')
    error('Site labels must be consecutive natural numbers starting from 1.');
end


thres = option.thres;
mite = option.mite;
K = option.K;
output = option.output;

disp(strcat('Edge screening rule: R-squared >', {' '}, string(thres), '%'))
disp(strcat('Number of MCMC iterations: ', {' '}, string(mite)))

if isnumeric(K)
    disp(strcat('Number of clusters specified as', {' '}, string(K)))
elseif strcmp(K, 'auto')
    disp('Number of clusters determined by the built-in method')
end

fprintf('Output path: %s\n', output);

%%
addpath('helpers'); 

% Construct X matrix from age and sex, and normalize FC and X. 
% Transfer site label into a matrix indicating batch indices.
preprocessed_data = process_data(FC, age, sex, site_label, output);

FC_norm = preprocessed_data.FC_norm;
X_norm = preprocessed_data.X_norm;
ind_batch = preprocessed_data.ind_batch;
data_info = preprocessed_data.data_info;

%%
% Gam0 is the matrix indicating the prescreen results,
% Gam0(i,j) = 1 indicates edge (i,j) has FC with R_squared above thres/100,
% thus included in the analysis
[Gam0, indepReg_info] = prescreen(FC_norm, X_norm, ind_batch, data_info, thres, output);

%%
% Parameter Initialization
% Specify hyperparameters
rho0 = 10^-2;
varrho0 = 10^-6;
initial_values = ...
    parameter_init(FC_norm, X_norm, ind_batch, Gam0, K, indepReg_info.FC_test_pred_mat, output);


%% MCMC
[MCMC_outputs, ~] = ...
    Gibbs_sampling(mite, X_norm, initial_values, ind_batch, rho0, varrho0, output);

end

