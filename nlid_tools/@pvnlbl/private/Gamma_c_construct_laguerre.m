function Gamma_c = Gamma_c_construct_laguerre(c,q,m,n,p)
Gamma_c = zeros(q*m*n*p,n*p);

cnt = 0;
ind2 = 0;
for i = 1:n
    for j = 1:p
        cnt = cnt+1;
        ind1 = ind2 + 1;
        ind2 = cnt*q*m;
        Gamma_c(ind1:ind2,cnt) = c;
    end
end
