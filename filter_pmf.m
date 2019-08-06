%Online PMF by (Salakhutdinov and Mnih, 2008)
% -dim: low rank
% -pen_u: penalty on u
% -pen_v: penalty on v

%Author: San Gultekin, Columbia University, san.gultekin@gmail.com



function PRED = filter_pmf(dim,pen_u,pen_v)
global DATA_MAT;
global ID_MAT;

T = size(DATA_MAT,2);
PRED = zeros(size(ID_MAT));

U = (2/dim)*rand(dim,size(DATA_MAT,1));
v = (2/dim)*rand(dim,1);

for t = 1:T
    %active measurements
    pl = ID_MAT(:,t);
    
    %forecasting
    if t>1
        PRED(:,t) = U(:,pl)'*v;
    end
    
    %current observation
    obs = DATA_MAT(pl,t);
    
    %Update (E-step)
    bar_U = U(:,pl);
    for ite=1:15
        sub_U = (pen_u*eye(dim) + v*v')\(pen_u*bar_U + v*obs');
        v = (pen_v*eye(dim) + sub_U*sub_U')\(pen_v*v + sub_U*obs);
    end  
    U(:,pl) = sub_U;
end
end