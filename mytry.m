clc
clear

date_=clock;
resultsfilename=['Results/r_Results_' num2str(date_(1)) '_' num2str(date_(2)) '_' num2str(date_(3)) '-' num2str(date_(4))  '_' num2str(date_(5)) '.txt'];
diary(resultsfilename);

mfilename('fullpath')                           % print filename of currently running script

load('DV_IJ_sparse_combine.mat')
load('X_m_neighbour_adoption_times_similarity_combine.mat')
y = DV_IJ_sparse_combine;
X = X_m_neighbour_adoption_times_similarity_combine;
clearvars DV_IJ_sparse_combine X_m_neighbour_adoption_times_similarity_combine

beta_0 = [-8.9801    0.1]

[b, hessian, grad, standard_error, covariance_matrix, t_stat, exit_flag, output] = band_runbi_ll_try(X, y, beta_0);

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

diary off 