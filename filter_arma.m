%ARMA predictor
% -nord: model order
% -reg: regularizer

%Author: San Gultekin, Columbia University, san.gultekin@gmail.com



function PRED = filter_arma(nord,reg)
global DATA_MAT;
global ID_MAT;

T = size(DATA_MAT,2);
PRED = zeros(size(ID_MAT));

lhs = reg*eye(nord);
rhs = zeros(nord,1);
weig = zeros(nord,1);

for t=1:T
    %forecasting
    if t==1
        PRED(:,t) = 0;
    elseif t>1 && t<=nord+1
        PRED(:,t) = func_predfill(DATA_MAT(:,t-1),ID_MAT(:,t-1),ID_MAT(:,t));
    elseif t>nord+1
        BLK = func_fill(DATA_MAT(:,t-1:-1:t-nord),ID_MAT(:,t-1:-1:t-nord));
        pred = BLK*weig;
        PRED(:,t) = pred(ID_MAT(:,t));
    end
    
    %AR-weight update
    if t>=nord+1
        BLK = func_fill(DATA_MAT(:,t-1:-1:t-nord),ID_MAT(:,t-1:-1:t-nord));
        lhs = lhs + BLK'*BLK;
        rhs = rhs + BLK'*func_fill(DATA_MAT(:,t),ID_MAT(:,t));
        weig = lhs\rhs;
    end
end
end