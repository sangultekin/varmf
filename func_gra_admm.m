%ADMM solver for filter_gra.

%Author: San Gultekin, Columbia University, san.gultekin@gmail.com



function [r, e] = func_gra_admm(L, x, lam1, rho, no_iter)
lam2 = 1/(1+rho);
dim2 = size(L,2);
e = zeros(size(x));
y = zeros(size(x));
for i=1:no_iter
    r = (1/rho)*(L'*L + lam1*eye(dim2)) \ (L'*(rho*(x-e)-y));
    a = x - L*r - e;    
    if a > lam2
        e = a - lam2;
    elseif a < -lam2
        e = a + lam2;
    else
        e = zeros(size(x));
    end
    y = y + rho*(L*r + e - x);
end
end