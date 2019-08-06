%Fill a matrix or vector.
%Helper function for algorithms that cannot handle sparsity.

%Author: San Gultekin, Columbia University, san.gultekin@gmail.com



function fX = func_fill(DATA_MAT,ID_MAT)
fX = zeros(size(DATA_MAT));
for i=1:size(DATA_MAT,2)
    %find where to fill
    pl = true(size(DATA_MAT(:,i)));
    pl(ID_MAT(:,i)) = false;
    %fill column
    fX(~pl,i) = DATA_MAT(~pl,i);
    fX(pl,i) = mean(fX(~pl,i)); 
end 
end