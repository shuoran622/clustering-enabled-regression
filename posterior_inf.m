function [clusters, mean_B, B_all_sorted, mean_sigma_sq, mean_delta_sq] = posterior_inf(M_all, B_all, sigma_sq_all, delta_sq_all)
    

    mite = size(M_all,3) - 1;
    ite_start = floor(mite/2) + 2;
    ite_end = mite + 1;
    mean_m = mean(M_all(:,:,ite_start:ite_end),3);
    mean_sigma_sq = mean(sigma_sq_all(:,:,ite_start:ite_end),3);
    mean_delta_sq = mean(delta_sq_all(:,:,:,ite_start:ite_end),4);

    mean_m_double = double(mean_m > 0.5); 
    
    %% Sort clusters by size
    disp('Sort region clusters by size and exclude empty clusters...')
    size_list = sum(mean_m_double,2);
    [size_list_sort,id_sort] = sort(size_list,'descend');

    K_list = id_sort(size_list_sort > 0); 

    fprintf('%d nonempty region clusters identified \n', length(K_list))
    B_all_sorted = B_all(K_list,K_list,:,:);
    mean_B = mean(B_all_sorted(:,:,:,ite_start:ite_end),4);
    
    mean_m_double = mean_m_double(K_list,:);
    
    
    clusters = {};
    for i = 1:length(K_list)
        clusters{i} = find(mean_m_double(i,:) == 1);
        fprintf('Cluster %d: ', i);
        disp(clusters{i});
    end

end
    