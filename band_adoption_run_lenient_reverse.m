% modified from band_adoption_run_ipc.m

clc
% clear
clear all

% disp('run complete cases with mean imputation and log(y_hat) as predicted trends and baseline probability as controls') 

load('matp_friend_reverse.mat');

% matp = [member_p friend_p band_p timeobs DV prob_adopt_week new_week_diff A_week_ijt];

X = matp(:,[6 7 8]);
X_N = N_prep_for_matp_jobs_organize;
S = combi_similarities;
% X = matp(:,[6 7]);
y = matp(:,5);

clearvars matp N_prep_for_matp_jobs_organize combi_similarities

% load('matpstd2.mat');                  % add old baseline probability to the model (two controls here: trend and BP)         mar 2017
% load('matp_bp_strict_rm.mat');
% load('newbp_store_rm_all_member_lenient.mat')       % variable name baseline_prob_store
% X = [X baseline_prob_store];
% clearvars baseline_prob_store

%         intercept gamma1     gamma2  baseline_prob   new_week_diff
beta_0 = [-5.0318    0.1763    2.0272    0.3235    0.0454];

[b, hessian, grad, standard_error, covariance_matrix, t_stat, exit_flag, output] = band_runbi_ll_reverse(X, y, beta_0);

save('b.mat','b') ;
save('standard_error.mat','standard_error') ;
save('t_stat.mat','t_stat') ;
save('exit_flag.mat','exit_flag') ;

display(b)
% display(standard_error)
display(t_stat)
display(grad)
display(output)

m = grad'*(-inv(hessian))*grad;          % convergence criterion of Train
