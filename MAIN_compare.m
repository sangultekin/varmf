%Code for Online Forecasting Matrix Factorization.
%Reference paper:
%S. Gultekin and J. Paisley. Nonlinear Kalman filtering with divergence minimization, 
%IEEE Transactions on Signal Processing, vol. 67, no. 5, pp. 1223-1236, 2019.

%Author: San Gultekin, Columbia University, san.gultekin@gmail.com



close all;clear all;clc; 
%#ok<*NUSED> 
%#ok<*NASGU>
home = pwd;

%load dataset
global DATA_MAT;  
global ID_MAT; 
cd('../data_parse');
load('DATA_MAT.mat');
cd(home);

%all sparsities
nnz_pct_arr = [10:10:100];

%parameter configuration
var_proc = 1e-4;
var_meas = 1e-4;
nord  = 24;
reg   = 1;
dim   = 5;
pen_u = 1;
pen_v = 1e-4;
tol_u = 1e-2;
%for rpca
lam1 = 1;lam2 = 1e-6;

%go
for nnz_ix=1:length(nnz_pct_arr)  
    %select sparse dataset here
    nnz_pct = nnz_pct_arr(nnz_ix);

    %get MC count
    cd('../data_parse');
    cd(['id_' num2str(nnz_pct)]);
    no_mc = length(dir)-2;
    cd(home);
    
    %run all algorithms
    for i=1:no_mc
        cd(['../data_parse/' 'id_' num2str(nnz_pct)]);
        load(['ID_MAT_' num2str(i)]);
        cd(home);
        
        PRED = filter_base();
        cd(['../results/compare_' num2str(nnz_pct)]);
        save(['PRED_base_' num2str(i) '.mat'],'PRED');
        cd(home);
        
        PRED = filter_arma(nord,reg);
        cd(['../results/compare_' num2str(nnz_pct)]);
        save(['PRED_arma_' num2str(i) '.mat'],'PRED');
        cd(home);
        
        PRED = filter_ckf(dim,var_proc,var_meas);
        cd(['../results/compare_' num2str(nnz_pct)]);
        save(['PRED_ckf_' num2str(i) '.mat'],'PRED');
        cd(home);
        
        PRED = filter_pmf(dim,pen_u,pen_v);
        cd(['../results/compare_' num2str(nnz_pct)]);
        save(['PRED_pmf_' num2str(i) '.mat'],'PRED');
        cd(home);
        
        PRED = filter_gra(dim,pen_u,pen_v);
        cd(['../results/compare_' num2str(nnz_pct)]);
        save(['PRED_gra_' num2str(i) '.mat'],'PRED');
        cd(home);

        PRED = filter_rpca(dim,pen_u,pen_v);
        cd(['../results/compare_' num2str(nnz_pct)]);
        save(['PRED_rpca_' num2str(i) '.mat'],'PRED');
        cd(home);
        
        PRED = filter_nmf(dim,pen_u,pen_v);
        cd(['../results/compare_' num2str(nnz_pct)]);
        save(['PRED_nmf_' num2str(i) '.mat'],'PRED');
        cd(home);
        
        PRED = filter_FP(dim,pen_u,pen_v,nord,reg);
        cd(['../results/compare_' num2str(nnz_pct)]);
        save(['PRED_FP_' num2str(i) '.mat'],'PRED');
        save(['STORE_FP_' num2str(i) '.mat'],'STORE');
        cd(home);
        
        PRED = filter_FT(dim,tol_u,pen_v,nord,reg);
        cd(['../results/compare_' num2str(nnz_pct)]);
        save(['PRED_FT_' num2str(i) '.mat'],'PRED');
        save(['STORE_FT_' num2str(i) '.mat'],'STORE');
        cd(home);
        
        PRED = filter_ZT(dim,pen_v,nord,reg);
        cd(['../results/compare_' num2str(nnz_pct)]);
        save(['PRED_LN_' num2str(i) '.mat'],'PRED');
        save(['STORE_LN_' num2str(i) '.mat'],'STORE');
        cd(home);
    
        disp(['Done processing: ' num2str(nnz_pct) ' | ' num2str(i)]);
    end
end