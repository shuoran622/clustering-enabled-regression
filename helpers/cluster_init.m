function [Initial_Cluster, optimalK] = cluster_init(FC_test_pred_mat, K)

% Calculate corr
d = size(FC_test_pred_mat, 1);
corr_mat = NaN(d,d);
for i = 1:d
    for j = 1:d
        
        if i <= j
            continue
        end
        
        corr_vec = [];
        for k = 1:d
            if k==i || k==j
                continue
            end         
            y1 = reshape(FC_test_pred_mat(max(i,k),min(i,k),:),[],1);
            y2 = reshape(FC_test_pred_mat(max(j,k),min(j,k),:),[],1);
            corr_vec = [corr_vec; corr(y1,y2)];
                
        end
        corr_mat(i,j) = mean(corr_vec);    
    end
end

if ischar(K)
    % NB customized
    tmp_corr_mat = corr_mat(~isnan(corr_mat));
    tmp_q = quantile(tmp_corr_mat, 0.95);
    A = double(corr_mat > tmp_q);
    A = A + A';
    D = diag(sum(A, 2));

    B = zeros(2*d,2*d);
    B((d+1):(2*d),1:d) = -eye(d);
    B(1:d,(d+1):(2*d)) = D-eye(d);
    B((d+1):(2*d),(d+1):(2*d)) = A;

    eigenvalues_B = eig(B);
    eigenvalues_B_real2 = NaN(2*d,1);
    for i = 1:(2*d)

        if abs(imag(eigenvalues_B(i))) < 10^-4
            eigenvalues_B_real2(i) = real(eigenvalues_B(i));
        else
            eigenvalues_B_real2(i) = -Inf;
        end

    end
    optimalK = sum(eigenvalues_B_real2 > 10^-4);
    fprintf('The initial optimal K is %d.\n', optimalK);
else
    optimalK = K;
    fprintf('The initial K is prespecified by the user as %d.\n', optimalK);
end
% spectral clustering
corr_mat_NonNeg = (corr_mat + 1) / 2;
corr_mat_NonNeg(isnan(corr_mat_NonNeg)) = 0;
W = corr_mat_NonNeg+corr_mat_NonNeg.';
W(W==0) = eps; % to ensure there's no degree issue
W = W - diag(diag(W)); % to ensure diagonal is 0

[dd,~] = size(W);

D = diag(sum(W, 2));

D_inv_sqrt = diag(1 ./ sqrt(diag(D)));
L_sym = eye(dd) - D_inv_sqrt * W * D_inv_sqrt;
L_sym = (L_sym + L_sym')/2;


[eigen_vector,~] = eigs(L_sym,d,'smallestabs'); 

eigen_vector = real(eigen_vector);

kmeans_idx = kmeans(normr(eigen_vector(:,1:optimalK)),optimalK,'Replicates',500);

Initial_Cluster = zeros(optimalK,d);
for i=1:d
    Initial_Cluster(kmeans_idx(i),i)=1;
end

end