function [R0, R] = Update_R(sigma_sq, delta_sq, Gam, ind_batch)

% R: S*d*d % reduced the size to decrease the computational burden 
% R0: S*d*d, R = sigma_sq.*R0

[S,num_batch] = size(ind_batch);
[d,~] = size(sigma_sq);

R0 = NaN(S,d,d);
R = NaN(S,d,d);
for i = 2:d
    for j = 1:(i-1)
        
        if Gam(i,j) == 1
            
            tmp_R0 = ind_batch(:,num_batch); 
            for num_b = 1:(num_batch-1)
                tmp_R0 = tmp_R0 + delta_sq(i,j,num_b).*ind_batch(:,num_b);
            end
                        
            tmp_R = sigma_sq(i,j).*tmp_R0;

            R0(:,i,j) = tmp_R0;
            R0(:,j,i) = tmp_R0;

            R(:,i,j) = tmp_R;
            R(:,j,i) = tmp_R;
                      
        end
        
    end
end


end