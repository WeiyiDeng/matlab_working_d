clc
% clear
clear all

% diary('wwdiary.txt')

date_=clock;
resultsfilename=['Results/r_Results_' num2str(date_(1)) '_' num2str(date_(2)) '_' num2str(date_(3)) '-' num2str(date_(4))  '_' num2str(date_(5)) '.txt'];
diary(resultsfilename);

% load('matp_b.mat');
% load('matp.mat');

% I = size(matp,1);
% numc = 3;                        % number of chunks-1
% lenc = fix(I/numc);       % length of each chunk

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
% X = matp(:,[6 7]);
% y = matp(:,5);
% X = zeros(lenc,26);

clearvars matp week1_dummy week2_dummy week3_dummy week4_dummy week5_dummy week6_dummy week7_dummy week8_dummy

% load('innov_contin.mat');
% innov_X = innov_contin;
% load('explor_contin.mat')

% dummy_X = [member_dummies member_dummies_week_d];
% X = mat1(:,[5 7]);
% y = mat1(:,4);
% beta_0 = zeros(1,size(X,2)+size(dummy_X,2)+1);
% beta_0 = zeros(1,size(X,2)+1);
% beta_0 = zeros(1,size(X,2)+2);
% beta_0 = [0 1 2 0 -0.007]
% beta_0 = [-6.1646   1    27.4751    3.0197   12.7780  zeros(1, size(innov_X,2)*2)]
beta_0 = [-6.1646   1    27.4751    3.0197   12.7780  zeros(1, 24)]
% beta_0 = [-6.1646   1    27.4751    3.0197   12.7780]
% beta_0 = [-6.1668   27.4752    3.0187   12.7781]
% beta_0 = [-100 100 -10]
% beta_0 = [-6.1474    3.1608    0.4154]
% [b, hessian, grad, standard_error, covariance_matrix, t_stat, exit_flag, output] = band_runbi_ll_p(X, y, dummy_X, beta_0);
% [b, hessian, grad, standard_error, covariance_matrix, t_stat, exit_flag, output] = band_runbi_ll_i(X, y, beta_0);
[b, hessian, grad, standard_error, covariance_matrix, t_stat, exit_flag, output] = band_runbi_ll_i(beta_0);

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
   
    