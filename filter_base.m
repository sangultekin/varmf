%Baseline predictor
% -observed entry: predict past value
% -missing entry: predict current average

%Author: San Gultekin, Columbia University, san.gultekin@gmail.com



function PRED = filter_base()
global DATA_MAT;
global ID_MAT;

T = size(DATA_MAT,2);
PRED = zeros(size(ID_MAT));
for t=2:T
    PRED(:,t) = func_predfill(DATA_MAT(:,t-1),ID_MAT(:,t-1),ID_MAT(:,t));
end
end