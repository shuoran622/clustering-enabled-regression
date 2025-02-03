function [] = plot_traj(k1, k2, B_all_sorted, data_info)

    mu = data_info.mean_X;
    sigma = data_info.std_X;

    mite =  size(B_all_sorted,4) - 1;
    ite_start = floor(mite/2) + 2;
    ite_end = mite + 1;
    num_samples = (ite_end - ite_start + 1);

    q = size(B_all_sorted, 3);

    age_min = floor(min(data_info.X(:,1)));
    age_max = ceil(max(data_info.X(:,1)));

    n_test = 100;
    tmp_x = linspace(age_min, age_max, n_test); 

    y_F_all = [];
    y_M_all = [];


    b_all = B_all_sorted(k1,k2,:,ite_start:ite_end);
    b_all = reshape(b_all,[q,num_samples]);

    beta0_F_all = NaN(num_samples,1);
    beta0_M_all = NaN(num_samples,1);
    beta1_F_all = NaN(num_samples,1);
    beta1_M_all = NaN(num_samples,1);
    beta2_F_all = NaN(num_samples,1);
    beta2_M_all = NaN(num_samples,1);
    for ite = 1:num_samples

        b = b_all(:,ite);

        beta0_F_all(ite) = - (b(1)*mu(1)/sigma(1) + b(2)*mu(2)/sigma(2) + b(3)*mu(3)/sigma(3) + ...
                b(4)*mu(4)/sigma(4) + b(5)*mu(5)/sigma(5));
        beta0_M_all(ite) = b(3)/sigma(3) + beta0_F_all(ite);


        beta1_F_all(ite) = b(1)/sigma(1);
        beta1_M_all(ite) = b(4)/sigma(4) + beta1_F_all(ite);

        beta2_F_all(ite) = b(2)/sigma(2);
        beta2_M_all(ite) = b(5)/sigma(5) +  beta2_F_all(ite);

    end

    beta_F_all = [beta0_F_all beta1_F_all beta2_F_all];
    beta_M_all = [beta0_M_all beta1_M_all beta2_M_all];


    tmp_y_F_all = beta_F_all(:,1) + beta_F_all(:,2)*(tmp_x) + ...
            beta_F_all(:,3)*(tmp_x).^2;
    tmp_y_M_all = beta_M_all(:,1) + beta_M_all(:,2)*(tmp_x) + ...
            beta_M_all(:,3)*(tmp_x).^2;

    y_F_u_all = quantile(tmp_y_F_all,0.975,1);
    y_F_l_all = quantile(tmp_y_F_all,0.025,1);

    y_M_u_all = quantile(tmp_y_M_all,0.975,1);
    y_M_l_all = quantile(tmp_y_M_all,0.025,1);

    %% Trajectories Plot
    linewidth = 4;

    figure('Position',[0 0 720 595]);
    plot(NaN, NaN, 'Color', 'b',...
        'LineWidth', linewidth, 'LineStyle', ':')
    hold on;
    plot(NaN, NaN,  'Color', 'r',...
        'LineWidth', linewidth);
    xlabel('Age');
    ylabel('FC');  
    legend({'Male', 'Female'});

    title('95% CI')
    plot(tmp_x, y_M_u_all, 'Color', 'b',...
        'LineWidth', linewidth, 'LineStyle', ':')
    hold on;
    plot(tmp_x, y_F_u_all,  'Color', 'r',...
        'LineWidth', linewidth);
    hold on
    plot(tmp_x, y_M_l_all, 'Color', 'b',...
        'LineWidth', linewidth, 'LineStyle', ':')
    hold on;
    plot(tmp_x, y_F_l_all,  'Color', 'r',...
        'LineWidth', linewidth);

    legend({'Male', 'Female'}, 'Location', 'best');
    legend boxoff
    
    set(gca, 'FontSize',  35)

end
