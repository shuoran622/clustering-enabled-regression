function sigma_sq = Update_sigma_sq(FC,X,B,M,R0,Gam,varrho0)

% permute FC to S*d*d
FC = permute(FC,[3 1 2]);

[S,d,~] = size(FC);
[~,q] = size(X);


% tmp_b
tmp_b_summary = NaN(q,d,d);
for l = 1:q
    tmp_b_summary(l,:,:) = M'*B(:,:,l)*M; %K*d
end

sigma_sq = NaN(d,d);

a = S/2 + varrho0;

for i = 2:d
    
    for j = 1:(i-1)
        
        if Gam(i,j) == 1
            tmp = FC(:,i,j) - X*tmp_b_summary(:,i,j);
        
            b = 1/2.*(tmp./R0(:,i,j)).'*tmp + varrho0;

            sigma_sq(i,j) = 1 / gamrnd(a, 1/b);
        end

    end
    
end


end