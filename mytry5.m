clc
clear

% date_=clock;
% resultsfilename=['Results/r_Results_' num2str(date_(1)) '_' num2str(date_(2)) '_' num2str(date_(3)) '-' num2str(date_(4))  '_' num2str(date_(5)) '.txt'];
% diary(resultsfilename);

% diary('wdiary.txt')

% mfilename('fullpath')                           % print filename of currently running script

load('DV_IJ_sparse_combine.mat')
load('X_m_neighbour_adoption_times_similarity_combine.mat')
load('X_baselineProbListens.mat');
load('X_baselineProbAdopts.mat');
load('X_predict_trend_log4061_original.mat')

y = DV_IJ_sparse_combine;
X = [X_m_neighbour_adoption_times_similarity_combine X_baselineProbListens X_baselineProbAdopts X_predict_trend_log4061_original];
% y = DV_IJ_sparse_combine(1:100000);
% X = [X_m_neighbour_adoption_times_similarity_combine(1:100000) X_baselineProbListens(1:100000) X_baselineProbAdopts(1:100000) X_predict_trend_log4061_original(1:100000)];
clearvars DV_IJ_sparse_combine X_m_neighbour_adoption_times_similarity_combine X_baselineProbListens X_baselineProbAdopts X_predict_trend_log4061_original

% X(100000:end,:) = [];
% y(100000:end) = [];

clearvars -EXCEPT X y

load('X_m_friend_adoption_times_similarity_combine.mat');
load('X_m_friend_adoption_times_similarity_combineD1.mat');
load('X_m_friend_adoption_times_similarity_combineD4.mat');

% X = [X X_m_friend_adoption_times_similarity_combine X_m_friend_adoption_times_similarity_combineD1 X_m_friend_adoption_times_similarity_combineD4];
% X = [X X_m_friend_adoption_times_similarity_combine(1:100000) X_m_friend_adoption_times_similarity_combineD1(1:100000) X_m_friend_adoption_times_similarity_combineD4(1:100000)];
X = [X X_m_friend_adoption_times_similarity_combine X_m_friend_adoption_times_similarity_combineD4];
% X = [X X_m_friend_adoption_times_similarity_combine(1:100000) X_m_friend_adoption_times_similarity_combineD4(1:100000)];

clearvars X_m_friend_adoption_times_similarity_combine X_m_friend_adoption_times_similarity_combineD1 X_m_friend_adoption_times_similarity_combineD4

% beta_0 = [-8.9801    0.1   0.1  0.1  0.1]
% beta_0 = [-9.2778    0.1021   -0.0004    0.0538    0.3893     0.1    0.1    0.1]
beta_0 = [-9.2778    0.1021   -0.0004    0.0538    0.3893     0.1     0.1]

disp(size(X))

[b, hessian, grad, standard_error, covariance_matrix, t_stat, exit_flag, output] = band_runbi_ll_try3(X, y, beta_0);

save('b.mat','b')
save('standard_error.mat','standard_error')
save('t_stat.mat','t_stat')
save('exit_flag.mat','exit_flag')

display(b)
% display(standard_error)
display(t_stat)
display(grad)
display(output)

m = grad'*(-inv(hessian))*grad;          % convergence criterion of Train

% diary off 