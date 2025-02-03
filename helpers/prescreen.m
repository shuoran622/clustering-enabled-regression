function [Gam0, indepReg_info] = prescreen(FC_norm, X_norm, ind_batch, data_info, thres)
    
    [d,~,S] = size(FC_norm);
    
    age_min = floor(min(data_info.X(:,1)));
    age_max = ceil(max(data_info.X(:,1)));
    
    n_test = 100;
    age_test = linspace(age_min, age_max, n_test).'; 
    X_test_female = [age_test, age_test.^2, zeros(n_test,1), zeros(n_test,1), zeros(n_test,1)];
    X_test_male = [age_test, age_test.^2, ones(n_test,1), age_test, age_test.^2];
    X_test = [X_test_female; X_test_male];
    X_test_norm = (X_test - data_info.mean_X)./(data_info.std_X);

    num_sites = size(ind_batch,2);
    id_sites = cell(num_sites,1);
    for nr = 1:num_sites
        id_sites{nr} = find(ind_batch(:,nr) == 1);
    end

    FC_permu = permute(FC_norm,[3 1 2]); % easy-to-use
    Gam0 = NaN(d,d); 
    R_squared = NaN(d,d);
    F_pvalue = NaN(d,d);
    FC_test_pred_mat = NaN(d,d,size(X_test_norm, 1));

    for i = 1:d
        for j = 1:d

            if i <= j
                continue
            end

            tmp_FC_ij = FC_permu(:,i,j);

            % OLS first to get weights
            tmp_md = fitlm(X_norm,tmp_FC_ij,'Intercept',false);
            tmp_e = tmp_md.Residuals(:,1);
            tmp_e = table2array(tmp_e);

            w_sites = NaN(num_sites,1);
            for nr = 1:num_sites
                w_sites(nr) = 1/nanmean(tmp_e(id_sites{nr}).^2);
            end

            w = NaN(S,1); 
            for ns = 1:S
                w(ns) = w_sites(ind_batch(ns,:) == 1);
            end

            % WLS
            tmp_md_w = fitlm(X_norm,tmp_FC_ij,'Intercept',false,'Weights',w);

            FC_test_pred = predict(tmp_md_w, X_test_norm);
            FC_test_pred_mat(i,j,:) = FC_test_pred;

            tmp_R_sq = tmp_md_w.Rsquared.Ordinary;

            if tmp_R_sq > thres/100
                Gam0(i,j) = 1;
            else
                Gam0(i,j) = 0;
            end

            R_squared(i,j) = tmp_R_sq;
            [F_pvalue(i,j),~] = coefTest(tmp_md_w);
        end
    end
    
    indepReg_info = struct('R_squared', R_squared, 'F_pvalue', F_pvalue, 'X_test', X_test, 'FC_test_pred_mat', FC_test_pred_mat);

end