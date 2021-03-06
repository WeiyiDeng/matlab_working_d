% modified from band_adoption_run_ipc.m

clc
% clear
clear all

global model_num newBP lenient

% diary('wwdiary.txt')

date_=clock;
resultsfilename=['Results/r_Results_' num2str(date_(1)) '_' num2str(date_(2)) '_' num2str(date_(3)) '-' num2str(date_(4))  '_' num2str(date_(5)) '.txt'];
diary(resultsfilename);

mfilename('fullpath')                           % print filename of currently running script


model_num = 1
newBP = 1
lenient = 1

if model_num == 1 && newBP == 0 && lenient == 0
    disp('run strict adoption with bands with no trend data deleted and log(y_hat) as predicted trends and baseline probability both as controls')
    load('matp_trend_strict_rm.mat');
    load('innov_contin_trend_strict_rm.mat');
    load('explor_contin_trend_strict_rm.mat');
elseif model_num == 1 && newBP == 1 && lenient == 0
    disp('run strict adoption with bands with no trend data deleted and log(y_hat) as predicted trends and new baseline probability both as controls')
    load('matp_trend_strict_rm.mat');
    load('innov_contin_trend_strict_rm.mat');
    load('explor_contin_trend_strict_rm.mat');
elseif model_num == 1 && newBP == 0 && lenient == 1
    disp('run lenient adoption with bands with no trend data deleted and log(y_hat) as predicted trends and baseline probability both as controls')
    load('matp_trend_lenient_rm.mat');
    load('innov_contin_trend_lenient_rm.mat');
    load('explor_contin_trend_lenient_rm.mat');
elseif model_num == 1 && newBP == 1 && lenient == 1
    disp('run lenient adoption with bands with no trend data deleted and log(y_hat) as predicted trends and new baseline probability both as controls')
    load('matp_trend_lenient_rm.mat');
    load('innov_contin_trend_lenient_rm.mat');
    load('explor_contin_trend_lenient_rm.mat');
elseif model_num == 2 && newBP == 0 && lenient == 0
    disp('Strict adoption no trend with baseline prob')
    load('matp_strict_adopt_std.mat');
    load('innov_contin2_strict_std.mat');
    load('explor_contin2_strict_std.mat');
elseif model_num == 2 && newBP == 1 && lenient == 0
    disp('Strict adoption no trend with new baseline prob with ad hoc smoothing')
    load('matp_strict_adopt_newBP_std.mat');
    load('innov_contin2_strict_std.mat');
    load('explor_contin2_strict_std.mat');
elseif model_num == 2 && newBP == 0 && lenient == 1
    disp('lenient adoption no trend with baseline prob')
    load('matp_lenient_adopt_std.mat');
    load('innov_contin2_lenient_std.mat');
    load('explor_contin2_lenient_std.mat');
elseif model_num == 2 && newBP == 1 && lenient == 1
    disp('lenient adoption no trend with new baseline prob with ad hoc smoothing')
    load('matp_lenient_adopt_newBP_std.mat');
    load('innov_contin2_lenient_std.mat');
    load('explor_contin2_lenient_std.mat');
elseif model_num == 3 && lenient == 0
    disp('run strict adoption with bands with no trend data deleted and log(y_hat) as predicted trends')
    load('matp_trend_strict_rm.mat');
    load('innov_contin_trend_strict_rm.mat');
    load('explor_contin_trend_strict_rm.mat');
elseif model_num == 3 && lenient == 1
    disp('run lenient adoption with bands with no trend data deleted and log(y_hat) as predicted trends')
    load('matp_trend_lenient_rm.mat');
    load('innov_contin_trend_lenient_rm.mat');
    load('explor_contin_trend_lenient_rm.mat');
end
    

% disp('run complete cases with mean imputation and exp(log(y_hat)) as predicted trends')                % print message
% disp('run complete cases with mean imputation and log(y_hat) as predicted trends') 
% disp('run complete cases with mean imputation and log(y_hat) as predicted trends and baseline probability as controls') 
% disp('Strict adoption no trend')
% disp('run strict adoption with bands with no trend data deleted and log(y_hat) as predicted trends')         % rm bands with no trend data
% disp('run strict adoption with bands with no trend data deleted and log(y_hat) as predicted trends and baseline probability both as controls') 

% % load('matp_b.mat');
% % load('matpstd2.mat');
% % load('innov_contin_trend.mat');
% % load('explor_contin_trend.mat');
% load('matp_trend.mat');                  % goes with the model with 30 variables (two trends)
% % load('matp_trend_exp.mat');
% load('innov_contin_std2.mat');
% load('explor_contin_std2.mat');
% load('matp_strict_adopt_std.mat');
% load('innov_contin2_strict_std.mat');
% load('explor_contin2_strict_std.mat');
% load('matp_trend_strict_rm.mat');
% load('innov_contin_trend_strict_rm.mat');
% load('explor_contin_trend_strict_rm.mat');

% week1_dummy = matp(:,7);
% week1_dummy(week1_dummy~=1)=0;
% 
% week2_dummy = matp(:,7);
% week2_dummy(week2_dummy==1)=0;
% week2_dummy(week2_dummy==2)=1;
% week2_dummy(week2_dummy~=1)=0;
% 
% week3_dummy = matp(:,7);
% week3_dummy(week3_dummy==1)=0;
% week3_dummy(week3_dummy==3)=1;
% week3_dummy(week3_dummy~=1)=0;
% 
% week4_dummy = matp(:,7);
% week4_dummy(week4_dummy==1)=0;
% week4_dummy(week4_dummy==4)=1;
% week4_dummy(week4_dummy~=1)=0;
% 
% week5_dummy = matp(:,7);
% week5_dummy(week5_dummy==1)=0;
% week5_dummy(week5_dummy==5)=1;
% week5_dummy(week5_dummy~=1)=0;
% 
% week6_dummy = matp(:,7);
% week6_dummy(week6_dummy==1)=0;
% week6_dummy(week6_dummy==6)=1;
% week6_dummy(week6_dummy~=1)=0;
% 
% week7_dummy = matp(:,7);
% week7_dummy(week7_dummy==1)=0;
% week7_dummy(week7_dummy==7)=1;
% week7_dummy(week7_dummy~=1)=0;
% 
% week8_dummy = matp(:,7);
% week8_dummy(week8_dummy==1)=0;
% week8_dummy(week8_dummy==8)=1;
% week8_dummy(week8_dummy~=1)=0;
% 
% week9_dummy = matp(:,7);
% week9_dummy(week9_dummy==1)=0;
% week9_dummy(week9_dummy==9)=1;
% week9_dummy(week9_dummy~=1)=0;
% 
% week10_dummy = matp(:,7);
% week10_dummy(week10_dummy==1)=0;
% week10_dummy(week10_dummy==10)=1;
% week10_dummy(week10_dummy~=1)=0;

% week1_dummy = matp(:,7);
% week1_dummy(week1_dummy==2)=1;
% week1_dummy(week1_dummy==3)=1;
% week1_dummy(week1_dummy==4)=1;
% week1_dummy(week1_dummy==5)=1;
% week1_dummy(week1_dummy~=1)=0;
% 
% week2_dummy = matp(:,7);
% week2_dummy(week2_dummy==1)=0;
% week2_dummy(week2_dummy==6)=1;
% week2_dummy(week2_dummy==7)=1;
% week2_dummy(week2_dummy==8)=1;
% week2_dummy(week2_dummy==9)=1;
% week2_dummy(week2_dummy==10)=1;
% week2_dummy(week2_dummy~=1)=0;
% 
% week3_dummy = matp(:,7);
% week3_dummy(week3_dummy==1)=0;
% week3_dummy(week3_dummy==11)=1;
% week3_dummy(week3_dummy==12)=1;
% week3_dummy(week3_dummy==13)=1;
% week3_dummy(week3_dummy==14)=1;
% week3_dummy(week3_dummy==15)=1;
% week3_dummy(week3_dummy~=1)=0;

% week4_dummy = matp(:,7);
% week4_dummy(week4_dummy==1)=0;
% week4_dummy(week4_dummy==4)=1;
% week4_dummy(week4_dummy==5)=1;
% week4_dummy(week4_dummy==6)=1;
% week4_dummy(week4_dummy==7)=1;
% week4_dummy(week4_dummy==8)=1;
% week4_dummy(week4_dummy~=1)=0;
% 
% week5_dummy = matp(:,7);
% week5_dummy(week5_dummy==1)=0;
% week5_dummy(week5_dummy==9)=1;
% week5_dummy(week5_dummy==10)=1;
% week5_dummy(week5_dummy==11)=1;
% week5_dummy(week5_dummy==12)=1;
% week5_dummy(week5_dummy==13)=1;
% week5_dummy(week5_dummy~=1)=0;
% 
% week6_dummy = matp(:,7);
% week6_dummy(week6_dummy==1)=0;
% week6_dummy(week6_dummy==14)=1;
% week6_dummy(week6_dummy==15)=1;
% week6_dummy(week6_dummy==16)=1;
% week6_dummy(week6_dummy==17)=1;
% week6_dummy(week6_dummy==18)=1;
% week6_dummy(week6_dummy~=1)=0;
% 
% week7_dummy = matp(:,7);
% week7_dummy(week7_dummy==1)=0;
% week7_dummy(week7_dummy==19)=1;
% week7_dummy(week7_dummy==20)=1;
% week7_dummy(week7_dummy==21)=1;
% week7_dummy(week7_dummy==22)=1;
% week7_dummy(week7_dummy==23)=1;
% week7_dummy(week7_dummy~=1)=0;
% 
% week8_dummy = matp(:,7);
% week8_dummy(week8_dummy==1)=0;
% week8_dummy(week8_dummy==24)=1;
% week8_dummy(week8_dummy==25)=1;
% week8_dummy(week8_dummy==26)=1;
% week8_dummy(week8_dummy==27)=1;
% week8_dummy(week8_dummy==28)=1;
% week8_dummy(week8_dummy~=1)=0;
% 
% week9_dummy = matp(:,7);
% week9_dummy(week9_dummy==1)=0;
% week9_dummy(week9_dummy==29)=1;
% week9_dummy(week9_dummy==30)=1;
% week9_dummy(week9_dummy==31)=1;
% week9_dummy(week9_dummy==32)=1;
% week9_dummy(week9_dummy==33)=1;
% week9_dummy(week9_dummy~=1)=0;
% 
% week10_dummy = matp(:,7);
% week10_dummy(week10_dummy==1)=0;
% week10_dummy(week10_dummy==34)=1;
% week10_dummy(week10_dummy==35)=1;
% week10_dummy(week10_dummy==36)=1;
% week10_dummy(week10_dummy==37)=1;
% week10_dummy(week10_dummy==38)=1;
% week10_dummy(week10_dummy~=1)=0;

% matp(:,7) = week1_dummy;

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
X = matp(:,[6 7]);
y = matp(:,5);

clearvars matp

if model_num == 1 && newBP == 0 && lenient == 0
    % load('matpstd2.mat');                  % add old baseline probability to the model (two controls here: trend and BP)         mar 2017
    load('matp_bp_strict_rm.mat');
    X = [X matp(:,6)];
    clearvars matp
elseif model_num == 1 && newBP == 0 && lenient == 1
    % load('matpstd2.mat');                  % add old baseline probability to the model (two controls here: trend and BP)         mar 2017
    load('matp_bp_lenient_rm.mat');
    X = [X matp(:,6)];
    clearvars matp
elseif model_num == 1 && newBP == 1 && lenient == 0
    % load('matpstd2.mat');                  % add new baseline probability to the model (two controls here: trend and BP)         mar 2017
    load('matp_newbp_strict_rm.mat');
    X = [X matp(:,6)];
    clearvars matp
elseif model_num == 1 && newBP == 1 && lenient == 1
    % load('matpstd2.mat');                  % add new baseline probability to the model (two controls here: trend and BP)         mar 2017
    load('matp_newbp_lenient_rm.mat');
    X = [X matp(:,6)];
    clearvars matp
end

% innov_X = innov_contin(:,1:4);
% % innov_X = [];
% explor_X = explor_contin(:,1:4);
% % explor_X = [];

innov_X = innov_contin;
explor_X = explor_contin;

clearvars innov_contin explor_contin
% dummy_X = [member_dummies member_dummies_week_d];
% X = mat1(:,[5 7]);
% y = mat1(:,4);
% beta_0 = zeros(1,size(X,2)+size(dummy_X,2)+1);
% beta_0 = zeros(1,size(X,2)+1);
% beta_0 = zeros(1,size(X,2)+2);
% beta_0 = [0 1 2 0 -0.007]
% beta_0 = [-6.1646   1    27.4751    3.0197   12.7780]
% beta_0 = [-6.1668   0.9337   27.4738    3.0426   12.7745    ones(1,2).*(-5)]
beta_0 = [-5.0318    0.1763    2.0272    0.3235    0.0454];
beta_0 = [beta_0         -0.0564    0.0017   -0.1336    0.0423   -0.0423   -0.0117    -0.0752    0.0117];  
beta_0 = [beta_0          0 0 0 0 0 0 0 0];
beta_0 = [beta_0         0 0 0 0 0 0 0 0];

if model_num == 1
    beta_0 = [beta_0 0]                                  % add one more beta for BP
else
end
% beta_0 = [beta_0         0.1    0.12   0    0.1   0.1   -0.1    0.1    -0.1];  
% beta_0 = [beta_0          0.04 -0.04 0.03 -0.03 -0.02 0.02 -0.03 0.01];
% beta_0 = [beta_0         -0.04 -0.01 0.01 0.01 0.03 0.01 0.01 -0.01]
% beta_0 = [beta_0         1    -1   -1    1   1   -1    1    -1];  
% beta_0 = [beta_0          1  1  1  1  1  1  1  1];
% beta_0 = [beta_0         1  1  1  1  1  1  1  1]
% beta_0 = [-5.0605   -3.5673   -1.8994    0.3361   -0.6996    0.0018    0.0510   -0.1177    0.0510   -0.0094   -0.0543    -0.0345    0.0065    0.9001    0.1132    0.8674    0.0227    0.8714    0.1070    0.9212    0.3865   -0.0094    -0.0047    0.0073   -0.0006    0.8617    0.1469    0.9941    0.5372]
% beta_0 = [-100 100 -10]
% beta_0 = [-6.1474    3.1608    0.4154]
% [b, hessian, grad, standard_error, covariance_matrix, t_stat, exit_flag, output] = band_runbi_ll_p(X, y, dummy_X, beta_0);
[b, hessian, grad, standard_error, covariance_matrix, t_stat, exit_flag, output] = band_runbi_ll_trend_format(X, y, beta_0, innov_X, explor_X);

display(b)
% display(standard_error)
display(t_stat)
display(grad)
display(output)

m = grad'*(-inv(hessian))*grad;          % convergence criterion of Train

betaCoefficient = b';
StandardError = standard_error;
tStatistics = t_stat';

if model_num == 1 && newBP == 0
    var_names = {'constant';'gamma_k';'gamma_*';'trend';'SI';'Innov_m';'Innov_m2';'Innov_f';'Innov_f2';...
        'Explor_m';'Explor_m2';'Explor_f';'Explor_f2';'Innov_m*Innov_f';'(Innov_m*Innov_f)2';...         % cell object here !
        'Explor_m*Explor_f';'(Explor_m*Explor_f)2';'Innov_m*SI';'Innov_m2*SI';'Innov_f*SI';'Innov_f2*SI';...
        'Explor_m*SI';'Explor_m2*SI';'Explor_f*SI';'Explor_f2*SI';'Innov_m*Innov_f*SI';'(Innov_m*Innov_f)2*SI';...
        'Explor_m*Explor_f*SI';'(Explor_m*Explor_f)2*SI';'baseline_prob'};
elseif model_num == 1 && newBP == 1
    var_names = {'constant';'gamma_k';'gamma_*';'trend';'SI';'Innov_m';'Innov_m2';'Innov_f';'Innov_f2';...
        'Explor_m';'Explor_m2';'Explor_f';'Explor_f2';'Innov_m*Innov_f';'(Innov_m*Innov_f)2';...         % cell object here !
        'Explor_m*Explor_f';'(Explor_m*Explor_f)2';'Innov_m*SI';'Innov_m2*SI';'Innov_f*SI';'Innov_f2*SI';...
        'Explor_m*SI';'Explor_m2*SI';'Explor_f*SI';'Explor_f2*SI';'Innov_m*Innov_f*SI';'(Innov_m*Innov_f)2*SI';...
        'Explor_m*Explor_f*SI';'(Explor_m*Explor_f)2*SI';'new_baseline_prob'};
elseif model_num == 2 && newBP == 0
    var_names = {'constant';'gamma_k';'gamma_*';'baseline_prob';'SI';'Innov_m';'Innov_m2';'Innov_f';'Innov_f2';...
        'Explor_m';'Explor_m2';'Explor_f';'Explor_f2';'Innov_m*Innov_f';'(Innov_m*Innov_f)2';...         % cell object here !
        'Explor_m*Explor_f';'(Explor_m*Explor_f)2';'Innov_m*SI';'Innov_m2*SI';'Innov_f*SI';'Innov_f2*SI';...
        'Explor_m*SI';'Explor_m2*SI';'Explor_f*SI';'Explor_f2*SI';'Innov_m*Innov_f*SI';'(Innov_m*Innov_f)2*SI';...
        'Explor_m*Explor_f*SI';'(Explor_m*Explor_f)2*SI'};
elseif model_num == 2 && newBP == 1
    var_names = {'constant';'gamma_k';'gamma_*';'new_baseline_prob';'SI';'Innov_m';'Innov_m2';'Innov_f';'Innov_f2';...
        'Explor_m';'Explor_m2';'Explor_f';'Explor_f2';'Innov_m*Innov_f';'(Innov_m*Innov_f)2';...         % cell object here !
        'Explor_m*Explor_f';'(Explor_m*Explor_f)2';'Innov_m*SI';'Innov_m2*SI';'Innov_f*SI';'Innov_f2*SI';...
        'Explor_m*SI';'Explor_m2*SI';'Explor_f*SI';'Explor_f2*SI';'Innov_m*Innov_f*SI';'(Innov_m*Innov_f)2*SI';...
        'Explor_m*Explor_f*SI';'(Explor_m*Explor_f)2*SI'};
elseif model_num == 3
    var_names = {'constant';'gamma_k';'gamma_*';'trend';'SI';'Innov_m';'Innov_m2';'Innov_f';'Innov_f2';...
        'Explor_m';'Explor_m2';'Explor_f';'Explor_f2';'Innov_m*Innov_f';'(Innov_m*Innov_f)2';...         % cell object here !
        'Explor_m*Explor_f';'(Explor_m*Explor_f)2';'Innov_m*SI';'Innov_m2*SI';'Innov_f*SI';'Innov_f2*SI';...
        'Explor_m*SI';'Explor_m2*SI';'Explor_f*SI';'Explor_f2*SI';'Innov_m*Innov_f*SI';'(Innov_m*Innov_f)2*SI';...
        'Explor_m*Explor_f*SI';'(Explor_m*Explor_f)2*SI'};
end

myTable = table(betaCoefficient, StandardError, tStatistics, 'RowNames', var_names)

mystr = strcat(num2str(date_(1)),num2str(date_(2)),num2str(date_(3)),num2str(date_(4)),num2str(date_(5)));
disp(mystr)

b_name = strcat('Results/mat_files/b_',mystr,'.mat');
SE_name = strcat('Results/mat_files/standard_error_',mystr,'.mat');
t_name = strcat('Results/mat_files/t_stat_',mystr,'.mat');
flag_name = strcat('Results/mat_files/exit_flag_',mystr,'.mat');

% save('b.mat','b') ;
% save('standard_error.mat','standard_error') ;
% save('t_stat.mat','t_stat') ;
% save('exit_flag.mat','exit_flag') ;

save(b_name,'b');
save(SE_name,'standard_error');
save(t_name,'t_stat');
save(flag_name,'exit_flag');

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
   
    