%Predict and fill the next observation vector.
%Helper function for algorithms that cannot handle sparsity.

%Author: San Gultekin, Columbia University, san.gultekin@gmail.com



function pred = func_predfill(data_prev,id_prev,id_curr)
    %find where to fill
    pl = true(size(data_prev));
    pl(id_prev) = false;
    %construct long prediction
    pred = zeros(size(data_prev));
    pred(~pl) = data_prev(~pl);
    pred(pl) = mean(pred(~pl));
    %construct short prediction
    pred = pred(id_curr,1);
end