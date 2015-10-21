clc
% clear
clear all

% diary('wwdiary.txt')

% format long g

date_=clock;
resultsfilename=['Results/r_Results_' num2str(date_(1)) '_' num2str(date_(2)) '_' num2str(date_(3)) '-' num2str(date_(4))  '_' num2str(date_(5)) '.txt'];
diary(resultsfilename);

% load('matp_b.mat');
% load('matp.mat');
% load('innov_contin.mat');
% load('explor_contin.mat');
load('matpstd2.mat');

%%
% X = matp(:,[6 7]);
% y = matp(:,5);
% 
% IVs = X;
% choice_dv = [y 1-y];
% 
% week_IV = 100*gampdf(IVs(:,2),0.8343,27.4712);              % w: NOTICE new lines here for fixed gamma parameters
% week_IV(IVs(:,2)<1)=0; 
% 
% % clearvars matp
% % load('innov_contin_std.mat');
% % load('explor_contin_std.mat');
% % 
% % prev_innov_IV = innov_contin(:,1:4);
% % prev_explor_IV = explor_contin(:,1:4);
% % clearvars innov_contin explor_contin
% % 
% % % prev_bs = [0.3900    0.1015    -0.0285    0.0031   -0.2426    0.2146    0.1718   -0.1603   -0.0624    0.0459]';
% % prev_bs = [0.3900    0.1015    -0.0509    0.0281   -0.2486    0.2183    0.1682   -0.1598    -0.0583    0.0421]';
% % prev_FV = [IVs(:,1) week_IV prev_innov_IV prev_explor_IV]*prev_bs;
% prev_bs = [0.3900    0.1015]';
% prev_FV = [IVs(:,1) week_IV]*prev_bs;
% 
% save('prev_FV.mat','prev_FV','-v7.3');
% clearvars prev_innov_IV prev_explor_IV
% 
% display('prev_FV')

% % w: draw a small random sample from the obs 
% temp1=rand(size(matp,1),1)>.985;     
% matp=matp(  matp(:,5)==1 | (matp(:,5)==0 &  temp1 ), :);
% matp=matp(rand(size(matp,1),1)>.95, :);
% sum(matp(:,5)==1)

% num = 200
% rand_vector = rand(num,1);
% y_simu = zeros(num,1);
% y_simu(rand_vector>0.5) = 1;
% 
% matp = zeros(num,7);
% matp(:,5) = y_simu;

% mat1 = csvread('imat4474.csv');
% load('member_dummies.mat');
% load('member_dummies_week_d.mat');

% matp(:,6) = matp(:,6)-mean(matp(:,6));
% matp(:,7) = matp(:,7)-mean(matp(:,7));
% matp(:,6) = (matp(:,6)-mean(matp(:,6)))./std(matp(:,6));
% matp(:,7) = (matp(:,7)-mean(matp(:,7)))./std(matp(:,7));

% X = [matp(:,[6 7]) week2_dummy week3_dummy week4_dummy week5_dummy week6_dummy week7_dummy week8_dummy week9_dummy week10_dummy];
% X = matp(:,[6 7 8]);

% load('matpstd.mat');
% 
X = matp(:,[6 7]);
y = matp(:,5);
clearvars matp
% 
% week_IV = 100*gampdf(X(:,2),0.8343,27.4712);              % w: NOTICE new lines here for fixed gamma parameters
% week_IV(X(:,2)<1)=0;   

% load('innov_contin_std2.mat');
% load('explor_contin_std2.mat');
% innov_IV = innov_contin(:,1:4);
% % innov_X = [];
% explor_IV = explor_contin(:,1:4);
% explor_X = [];
% innov_IV = innov_contin(:,5:6);
% prev_innov_IV = innov_contin(:,5:6)*[-0.0112   -0.0003]';

% clearvars innov_contin explor_contain  
% 
% % innov_WD_multip = zeros(size(innov_IV));
% % for i = 1:size(innov_IV,2);
% %     innov_WD_multip(:,i) = innov_IV(:,i).*week_IV;
% % end
% innov_WD_multip = [];
% clearvars innov_IV

% load('explor_contin_std.mat');

% explor_IV = explor_contin(:,5:6);
% prev_explor_IV = explor_contin(:,5:6)*[-0.0245    0.0172]';
% clearvars explor_contin
% 
% explor_WD_multip = zeros(size(explor_IV));
% for j = 1:size(explor_IV,2);
%     explor_WD_multip(:,j) = explor_IV(:,j).*week_IV;
% end
% explor_WD_multip = [];

% clearvars explor_IV

% dummy_X = [member_dummies member_dummies_week_d];
% X = mat1(:,[5 7]);
% y = mat1(:,4);
% beta_0 = zeros(1,size(X,2)+size(dummy_X,2)+1);
% beta_0 = zeros(1,size(X,2)+1);
% beta_0 = zeros(1,size(X,2)+2);
% beta_0 = [0 1 2 0 -0.007]
% beta_0 = [-5.6204    0.8343   27.4712    0.3900    0.1015]
% beta_0 = [-5.0318    1.1928    2.0272    0.3235    0.0454]
% beta_0 = [0    1    0    0    0]
beta_0 = [-5.0277    1.4538    1.5231    0.3256    0.0366]
% beta_0 = [-5.0354    0.6087   3.31313818317199    0.3250    0.0577]        % log(prev_b(3)) here
% beta_0 = [-6.1668   0.9337   27.4738    3.0426   12.7745    ones(1,2).*(-5)]
% beta_0 = [-0.0217   -0.0317    -0.2376    0.1782    0.0564   -0.1002    -0.1336    0.1002]
% beta_0 = [-1.5009    -0.2503    0.5233    -0.4710]
% beta_0 = [0.0170 0.0170 0.0070 0.0070]
% beta_0 = [-100 100 -10]
% beta_0 = [-6.1474    3.1608    0.4154]
% [b, hessian, grad, standard_error, covariance_matrix, t_stat, exit_flag, output] = band_runbi_ll_p(X, y, dummy_X, beta_0);
[b, hessian, grad, standard_error, covariance_matrix, t_stat, exit_flag, output] = band_runbi_ll_c(X, y, beta_0);
% [b, hessian, grad, standard_error, covariance_matrix, t_stat, exit_flag, output] = band_runbi_ll_i0_0(X, y, beta_0, innov_WD_multip, explor_WD_multip, prev_innov_IV, prev_explor_IV);

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

%%
% load('row_mid.mat');
% load('row_num.mat');
% 
% I = 158
% 
% b_store_p = cell(I,1);
% se_store_p = cell(I,1);
% cov_store_p = cell(I,1);
% tstat_store_p = cell(I,1);
% % hess = cell(I,1);
% i_id = zeros(I,1);
% est_exitflag_p = zeros(I,1);
% 
% ind = 0;
% for i = 1:I
%     i_id(i) = row_mid(i);
%     nrows = row_num(i);
%     imat = matp(ind+1:ind+nrows,:);         % ind tbc
%     
%     X = imat(:,[6 7]);
%     y = imat(:,5);
%     beta_0 = zeros(1,size(X,2)+1);
%     [b, hessian, standard_error, covariance_matrix, t_stat, exit_flag] = band_runbi_ll_p(X, y, beta_0);
% 
%     b_store_p{i} = b;
%     se_store_p{i} = standard_error';
%     cov_store_p{i} = covariance_matrix;
%     tstat_store_p{i} = t_stat;
% %    hess{i} = hessian;
%     est_exitflag_p(i) = exit_flag;
%         
%     ind = ind+nrows; 
% %     if i >= 10
% %         break
% %     end    
% end
%    
% save('b_store_p.mat','b_store_p') ;
% save('se_store_p.mat','se_store_p') ;
% save('tstat_store_p.mat','tstat_store_p') ;
% save('est_exitflag_p.mat','est_exitflag_p') ;
   
    