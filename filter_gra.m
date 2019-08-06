%Fixed penalty (FP) matrix factorization.
% -dim: low rank
% -pen_u: penalty on u
% -pen_v: penalty on v
% -nord: model order
% -reg: M-step regularizer

%Author: San Gultekin, Columbia University, san.gultekin@gmail.com



function PRED = filter_gra(dim,pen_u,pen_v)
global DATA_MAT;
global ID_MAT;

T = size(DATA_MAT,2);
PRED = zeros(size(ID_MAT));

U = (2/dim)*rand(dim,size(DATA_MAT,1));
V = (2/dim)*rand(dim,1);

for t = 1:T
    %active measurements
    pl = ID_MAT(:,t);

    PRED(:,t) = U(:,pl)'*V;
    
    %current observation
    obs = DATA_MAT(pl,t);
    
    %Update (E-step)
    for ite=1:15
        V = func_gra_admm(U(:,pl)', obs, pen_v, 1, 5);
        U(:,pl) = (pen_u*eye(dim) + V*V')\(V*obs');
    end    
end
end