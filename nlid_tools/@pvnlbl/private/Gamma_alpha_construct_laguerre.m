function Gamma_alpha = Gamma_alpha_construct_laguerre(alpha,q,m,n,p)
Gamma_alpha = zeros(q*m*n*p,q*m);

cnt = 0;
ind2 = 0;
for i = 1:n
    for j = 1:p
        cnt = cnt+1;
        row_block = diag(alpha(cnt,1)*ones(q*m,1));
        ind1 = ind2+1;
        ind2 = cnt*q*m;
        Gamma_alpha(ind1:ind2,:) = row_block;
    end
end


