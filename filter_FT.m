%Fixed tolerance (FT) matrix factorization.
% -dim: low rank
% -tol_u: penalty on u
% -pen_v: penalty on v
% -nord: model order
% -reg: M-step regularizer

%Author: San Gultekin, Columbia University, san.gultekin@gmail.com



function PRED = filter_FT(dim,tol_u,pen_v,nord,reg)
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
        hV = zeros(size(V));
    else
        Vn = BLK*weig;
        PRED(:,t) = U(:,pl)'*Vn;
        hV = Vn;
    end
    
    %current observation
    obs = DATA_MAT(pl,t);
    
    hU = U(:,pl); %previous U values
    if t==1; hU = hU * 0; end;
    for ite=1:15
        %update v
        V = (pen_v*eye(dim) + U(:,pl)*U(:,pl)')\(pen_v*hV + U(:,pl)*obs);
        %construct polynomial
        a1 = sum( obs.*(hU'*V) , 1);
        a2 = norm(V)^2;
        a3 = norm(obs)^2;
        a4 = norm(hU'*V)^2;
        a5 = tol_u^2 - norm(obs)^2;
        pol = [ -(a3+a5)*a2^2 , -2*(a3+a5)*a2 , a4-2*a1-a5 ];
        %find roots
        rt = real( roots(pol) ); %can have +i.0000
        lagm = max(rt);
        pen_u = 1/lagm;
        %update U
        U(:,pl) = (pen_u*eye(dim) + V*V')\(pen_u*hU + V*obs');
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