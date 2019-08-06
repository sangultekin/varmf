%CKF by (Gultekin and Paisley, 2014)
% -dim: low rank
% -var_proc: process noise variance
% -var_meas: measurement noise variance

%Author: San Gultekin, Columbia University, san.gultekin@gmail.com



function PRED = filter_ckf(dim,var_proc,var_meas)
global DATA_MAT;
global ID_MAT;

T = size(DATA_MAT,2);
PRED = zeros(size(ID_MAT));

obj_u.mu  = (2/dim)*rand(dim,size(DATA_MAT,1));
obj_u.Sig = 1*eye(dim);
obj_v.mu  = (2/dim)*rand(dim,1);
obj_v.Sig = 1*eye(dim);

for t = 1:T 
    %active measurements
    pl = ID_MAT(:,t);
    
    %forecasting
    if t>1
        PRED(:,t) = obj_u.mu(:,pl)'*obj_v.mu;
    end
    
    %current observation
    obs = DATA_MAT(pl,t);
    
    %prior 
    prior_mu1  = obj_u.mu(:,pl);
    prior_Sig1 = obj_u.Sig + var_proc*eye(dim);
    prior_mu2  = obj_v.mu;
    prior_Sig2 = obj_v.Sig + var_proc*eye(dim);
    
    %posterior
    post_mu1  = prior_mu1;
    post_Sig1 = prior_Sig1;
    post_mu2  = prior_mu2;
    post_Sig2 = prior_Sig2;     
    for ite=1:15
        post_Sig1 = (prior_Sig1\eye(dim) + (post_mu2*post_mu2' + post_Sig2)/var_meas) \ eye(dim);
        post_mu1  = post_Sig1 * (prior_Sig1\prior_mu1 + post_mu2*obs'/var_meas);
        post_Sig2 = (prior_Sig2\eye(dim) + (post_mu1*post_mu1' + size(ID_MAT,1)*post_Sig1)/var_meas) \ eye(dim);
        post_mu2  = post_Sig2 * (prior_Sig2\prior_mu2 + post_mu1*obs/var_meas);
    end  
    
    %update observed part
    obj_u.mu(:,pl) = post_mu1;
    obj_u.Sig      = post_Sig1;
    obj_v.mu       = post_mu2;
    obj_v.Sig      = post_Sig2;   
end
end