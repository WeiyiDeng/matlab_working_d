% modified from band_adoption_run_ipc.m

clc
% clear
clear all

% disp('run complete cases with mean imputation and log(y_hat) as predicted trends and baseline probability as controls') 

load('matp_friend_reverse.mat');
predict_trend = csvread('predict_trend_log4061_lenient.csv',1,0);

index = find(ismember(matp(:,3),predict_trend(:,1)));

predict_trend_mat = sparse(predict_trend(:,1),predict_trend(:,2),predict_trend(:,3),...
    max(predict_trend(:,1)),max(predict_trend(:,2)));

matp = matp(index,:);

trend_hat = zeros(size(matp,1),1);
for r = 1:length(trend_hat)
    trend_hat(r) = predict_trend_mat(matp(r,3),matp(r,4));
end
% matp(:,6) = trend_hat;

band_topics = csvread('band_count_topics10.csv',1,0);
band_topics_mat = sparse(band_topics(:,1),1,band_topics(:,2),max(band_topics(:,1)),1);
topics_count = zeros(size(matp,1),1);
for r = 1:length(topics_count)
    topics_count(r) = band_topics_mat(matp(r,3));
end

band_birth = csvread('introdate3.csv',1,0);
band_birth_mat = sparse(band_birth(:,1),1,band_birth(:,2),max(band_birth(:,1)),1);
band_age = zeros(size(matp,1),1);
for r = 1:length(band_age)
    band_age(r) = matp(r,4)-band_birth_mat(matp(r,3));
end
sum(band_age<0)

hist(band_topics(:,2))
% hist(band_birth(:,2))

% matp = [member_p friend_p band_p timeobs DV prob_adopt_week new_week_diff A_week_ijt];

X = matp(:,[6 7 8]);
% X = matp(:,[6 7]);
y = matp(:,5);

clearvars matp 

load('N_prep_for_matp_jobs_organize.mat');
load('combi_similarities.mat')

N_prep_for_matp_jobs_organize = N_prep_for_matp_jobs_organize(index,:);

[row col val] = find(N_prep_for_matp_jobs_organize==1);
store_r = zeros(length(row),1);
for i = 1:length(row)-1
    if row(i+1) == row(i)+2
        store_r(i) = 1;
    end
end
row0 = row(store_r==1)+1;
col0 = col(store_r==1);
val0 = zeros(size(col0));

row_store = [row; row0];
col_store = [col; col0];
val_store = [val; val0];

% val_pdf = normpdf(val_store);
% tic
% N_prep_for_matp_jobs_organize=sparse(row_store,col_store,val_pdf);
% toc

None0s_X_N = [row_store col_store val_store];
X_N = N_prep_for_matp_jobs_organize;
S = combi_similarities;

clearvars N_prep_for_matp_jobs_organize combi_similarities

% load('matpstd2.mat');                  % add old baseline probability to the model (two controls here: trend and BP)         mar 2017
% load('matp_bp_strict_rm.mat');
% load('newbp_store_rm_all_member_lenient.mat')       % variable name baseline_prob_store
% X = [X baseline_prob_store];
% clearvars baseline_prob_store

%         intercept gamma1    gamma2  baseline_prob gtrend(5)  new_week_diff  band_age  topics_count  Age*new_week_diff     Diversity*new_week_diff(10)     Diversity*Age    Diversity*Age*new_week_diff     week_diff_coef_n  normal_sd   similarity_coef
% beta_0 = [-5.0318    0.1763    2.0272    0.3235       0.1      0.0454          0.1         0.1           0.1                     0.1                              0.1                       0.1                      0.1            1          1];
%         intercept gamma1    gamma2  baseline_prob gtrend(5)  new_week_diff  band_age  topics_count  week_diff_coef_n  normal_sd(10)   similarity_coef
% beta_0 = [-5.0318    0.1763    2.0272    0.3235       0.1      0.0454          0.1         0.1             0.1            1              1];
beta_0 = [-5.0956    0.1766    2.0276    0.0401    0.0671    0.0176   -0.0065      -0.2234    0.0153    1.0260    1.0120]

[b, hessian, grad, standard_error, covariance_matrix, t_stat, exit_flag, output] = band_runbi_ll_trend_reverse_N_age_diverse_try(X, trend_hat, band_age, topics_count, X_N, None0s_X_N, S, y, beta_0);

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
