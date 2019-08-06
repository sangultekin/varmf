%Zero tolerance (ZT) matrix factorization.
% -dim: low rank
% -pen_v: penalty on v
% -nord: model order
% -reg: M-step regularizer

%Author: San Gultekin, Columbia University, san.gultekin@gmail.com



function PRED = filter_ZT(dim,pen_v,nord,reg)
global DATA_MAT;
global ID_MAT;

T = size(DATA_MAT,2);
PRED = zeros(size(ID_MAT));

U = (2/dim)*rand(dim,size(DATA_MAT,1));
V = (2/dim)*rand(dim,1);

%Past values of v for forecasting.
BLK = zeros(dim,nord);

lhs = reg*eye(nord);
rhs = zeros(nord,1);
weig = zeros(nord,1);

for t = 1:T
    %active measurements
    pl = ID_MAT(:,t);
    
    %forecasting
    if t==1
        PRED(:,t) = 0;
        hV = zeros(size(V));
    elseif t>1 && t<=nord+1
        PRED(:,t) = U(:,pl)'*V;
        hV = V;
    elseif t>nord+1
        Vn = BLK*weig;
        PRED(:,t) = U(:,pl)'*Vn;
        hV = Vn;
    end
    
    %current observation
    obs = DATA_MAT(pl,t);
     
    %Update (E-step)
    hU = U(:,pl); %previous U value
    for ite=1:15
        %update v
        V = (pen_v*eye(dim) + U(:,pl)*U(:,pl)')\(pen_v*hV + U(:,pl)*obs);
        %update U
        lagm = (hU'*V - obs)'/norm(V)^2;
        vrep = repmat(V,1,size(ID_MAT,1));
        U(:,pl) = hU - bsxfun(@times,lagm,vrep);
    end 
    
    %Update (M-step)
    if t>=nord+1       
        lhs = lhs + BLK'*BLK;
        rhs = rhs + BLK'*V;
        weig = lhs\rhs;
    end
    
    %Update BLK.
    BLK = [V BLK(:,1:end-1)];    
end
end