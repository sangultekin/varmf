%Projection function for filter_rpca.

%Author: San Gultekin, Columbia University, san.gultekin@gmail.com



function [r, e] = func_rpca_proj(L, x, lam1, lam2, no_iter)
dim2 = size(L,2);
e = zeros(size(x));

for i=1:no_iter
    r = (L'*L + lam1*eye(dim2)) \ (L'*(x-e));
    a = x - L*r;    
    if a > lam2
        e = a - lam2;
    elseif a < -lam2
        e = a + lam2;
    else
        e = zeros(size(x));
    end    
end
end