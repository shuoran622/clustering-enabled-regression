function p = Update_p(M, alpha)

tmp = sum(M,2) + alpha;

tmp = tmp';
P_tmp = gamrnd(repmat(tmp,1,1),1,1,length(tmp));

P_tmp = P_tmp ./ repmat(sum(P_tmp,2),1,length(tmp));

p = P_tmp';


end