function preprocessed_data =...
    process_data(FC, age, sex, site_label, output)

disp('---------- Preprocess data ----------')

[~,~,S] = size(FC);

ind_batch = zeros(S,length(unique(site_label)));
for s = 1:S
    ind_batch(s, site_label(s)) = 1;
end

disp(strcat('Create the batch indicator for sites:', {' '}, strjoin(string(unique(site_label')))))


FC_norm = normalize(FC, 3);
mean_FC = mean(FC, 3);
std_FC = std(FC, 0, 3);

disp('Standardize FC: Finished')


X = NaN(S,5);
X(:,1) = age; % age
X(:,2) = age.^2; % age square 
X(:,3) = sex;
X(:,4) = sex .* age;
X(:,5) = sex .* age.^2;


X_norm = normalize(X);
mean_X = mean(X, 1);
std_X = std(X, 1);


disp('Create and standardize covariates based on age and sex: Finished')


data_info = struct(...
    'mean_FC', mean_FC, ...
    'std_FC', std_FC, ...
    'mean_X', mean_X, ...
    'std_X', std_X, ...
    'X', X);

preprocessed_data = struct(...
    'FC_norm', FC_norm,...
    'X_norm', X_norm, ...
    'ind_batch', ind_batch, ...
    'data_info', data_info);

if ~exist(output, 'dir')
    mkdir(output);
end


save_path = strcat(output,'/preprocessed_data.mat');
save(save_path, 'preprocessed_data', '-v7.3');


fprintf('Preprocessed dat saved as %s\n', save_path);


disp('---------- Done ----------')


end