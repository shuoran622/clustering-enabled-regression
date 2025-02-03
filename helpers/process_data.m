function [FC_norm, X_norm, ind_batch, data_info] = process_data(FC, age, sex, site_label)

    [d,~,S] = size(FC);
    
    ind_batch = zeros(S,length(unique(site_label)));
    for s = 1:S
        ind_batch(s, site_label(s)) = 1;
    end
   
    X = NaN(S,5);
    X(:,1) = age; % age
    X(:,2) = age.^2; % age square 
    X(:,3) = sex;
    X(:,4) = sex .* age;
    X(:,5) = sex .* age.^2;
    
    FC_norm = normalize(FC, 3);
    mean_FC = mean(FC, 3);
    std_FC = std(FC, 0, 3);

    X_norm = normalize(X);
    mean_X = mean(X, 1);
    std_X = std(X, 1);
    
    data_info = struct('mean_FC', mean_FC, 'std_FC', std_FC, 'mean_X', mean_X, 'std_X', std_X, 'X', X);
end