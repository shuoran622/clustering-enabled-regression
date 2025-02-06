function data = load_data(file_name, var_names)

data = load(file_name, var_names{:});

disp('---------- Load the data ----------')
[num_region,~,num_subject] = size(data.FC);

disp(strcat('Number of subjects:', {' '}, string(num_subject)))
disp(strcat('Number of brain regions:', {' '}, string(num_region)))

disp('---------- Done ----------')


end