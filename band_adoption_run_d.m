clc
% clear
clear all

% diary('wwdiary.txt')

date_=clock;
resultsfilename=['Results/r_Results_' num2str(date_(1)) '_' num2str(date_(2)) '_' num2str(date_(3)) '-' num2str(date_(4))  '_' num2str(date_(5)) '.txt'];
diary(resultsfilename);

% load('matp_b.mat');
load('matp.mat');

week1_dummy = matp(:,7);
week1_dummy(week1_dummy==2)=1;
week1_dummy(week1_dummy==3)=1;
week1_dummy(week1_dummy==4)=1;
week1_dummy(week1_dummy~=1)=0;

% week2_dummy = matp(:,7);
% week2_dummy(week2_dummy==1)=0;
% week2_dummy(week2_dummy==5)=1;
% week2_dummy(week2_dummy==6)=1;
% week2_dummy(week2_dummy==7)=1;
% week2_dummy(week2_dummy==8)=1;
% week2_dummy(week2_dummy~=1)=0;
% 
% week3_dummy = matp(:,7);
% week3_dummy(week3_dummy==1)=0;
% week3_dummy(week3_dummy==9)=1;
% week3_dummy(week3_dummy==10)=1;
% week3_dummy(week3_dummy==11)=1;
% week3_dummy(week3_dummy==12)=1;
% week3_dummy(week3_dummy~=1)=0;
% 
% week4_dummy = matp(:,7);
% week4_dummy(week4_dummy==1)=0;
% week4_dummy(week4_dummy==13)=1;
% week4_dummy(week4_dummy==14)=1;
% week4_dummy(week4_dummy==15)=1;
% week4_dummy(week4_dummy==16)=1;
% week4_dummy(week4_dummy~=1)=0;
% 
% week5_dummy = matp(:,7);
% week5_dummy(week5_dummy==1)=0;
% week5_dummy(week5_dummy==17)=1;
% week5_dummy(week5_dummy==18)=1;
% week5_dummy(week5_dummy==19)=1;
% week5_dummy(week5_dummy==20)=1;
% week5_dummy(week5_dummy~=1)=0;
% 
% week6_dummy = matp(:,7);
% week6_dummy(week6_dummy==1)=0;
% week6_dummy(week6_dummy==21)=1;
% week6_dummy(week6_dummy==22)=1;
% week6_dummy(week6_dummy==23)=1;
% week6_dummy(week6_dummy==24)=1;
% week6_dummy(week6_dummy~=1)=0;

% week2_dummy = matp(:,7);
% % sum(week1_dummy==1)
% week1_dummy(week1_dummy~=1)=0;
matp(:,7) = week1_dummy;
% week2_dummy(week2_dummy~=2)=0;
% week2_dummy(week2_dummy==2)=1;
% matp(:,8) = week2_dummy;

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

load('member_dummies_week_d.mat');

% X = [matp(:,[6 7])  week2_dummy week3_dummy week4_dummy week5_dummy week6_dummy];
% X = [matp(:,[6 7])  member_dummies_week_d];
X = matp(:,[6 7]);
y = matp(:,5);

clearvars matp week1_dummy week2_dummy week3_dummy week4_dummy week5_dummy week6_dummy
dummy_X = member_dummies_week_d;
% X = mat1(:,[5 7]);
% y = mat1(:,4);
beta_0 = zeros(1,size(X,2)+size(dummy_X,2)+1);
% beta_0 = zeros(1,size(X,2)+1);
% beta_0 = zeros(1,size(X,2)+2);
% beta_0 = [0 1 2 0 -0.007]
% beta_0 = [0 1 2 0]
% beta_0 = [-100 100 -10]
% beta_0 = [-5.6814 2.6098 -0.0040]

% if matlabpool('size') == 0
%     matlabpool local 2
% end

[b, hessian, grad, standard_error, covariance_matrix, t_stat, exit_flag, output] = band_runbi_ll_d(X, y, dummy_X, beta_0);
% [b, hessian, grad, standard_error, covariance_matrix, t_stat, exit_flag, output] = band_runbi_ll_d(X, y, beta_0);

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

% matlabpool CLOSE

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
   
%% cut the wide (sparse) matrix into chunks
A = rand(8009,25);
b = [1:25]';
D = zeros(size(A,1),1);
numc = 10;                        % number of chunks
lenc = fix(size(A,1)/numc);       % length of each chunk
for i = 1:numc
    D((i-1)*lenc+1:i*lenc) = A((i-1)*lenc+1:i*lenc,:)*b;
end
if rem(size(A,1),numc)>0
    D(i*lenc:end) = A(i*lenc:end,:)*b;
else
end
C = A*b;
sum(D~=C)

    
