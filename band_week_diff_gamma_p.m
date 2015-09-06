clc
% clear
clear all

% diary('wwdiary.txt')

% date_=clock;
% resultsfilename=['Results/r_Results_' num2str(date_(1)) '_' num2str(date_(2)) '_' num2str(date_(3)) '-' num2str(date_(4))  '_' num2str(date_(5)) '.txt'];
% diary(resultsfilename);

% load('matp_b.mat');
load('matp.mat');

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

X = matp(:,[6 7 8]);
y = matp(:,5);

IVs = X;
choice_dv = [y 1-y];

% dummy_X = [member_dummies member_dummies_week_d];
% X = mat1(:,[5 7]);
% y = mat1(:,4);
% beta_0 = zeros(1,size(X,2)+size(dummy_X,2)+1);
% beta_0 = zeros(1,size(X,2)+1);
% beta_0 = zeros(1,size(X,2)+2);

b_k = 0:0.5:5;
b_theta = 0:1:10;
b_wd = 0:0.2:2;                   % 0.2495

display('LL')

ll_cube = zeros(length(b_k),length(b_theta),length(b_wd));
for i = 1:length(b_k)
    for j = 1:length(b_theta)
        for d = 1:length(b_wd)
            
            b = [-6.1474 b_k(i) b_theta(j) 3.1975 b_wd(d)];
            
            const = b(1);
            bs = b(4:end)';
            week_IV = IVs(:,3).*gampdf(IVs(:,2),b_k(i),b_theta(j));
            FV = [IVs(:,1) week_IV]*bs;
            
            exp_util = exp(-(const+FV));         % this is now the utility of the external good
            prob=1./(1+exp_util);                % this is still the probability of choosing the product
            pmat = [prob 1-prob];
            pmat = pmat.*choice_dv;
            [r c p] = find(pmat);                                             % I*1
            ll_cube(i,j,d) = sum(log(p));
        end
    end
end

mat_3d_name = ['ll_cube_' num2str(b_k(1)) '_' num2str(b_k(end)) '_' num2str(b_theta(1)) '_' num2str(b_theta(end)) '_' num2str(b_wd(1)) '_' num2str(b_wd(end)) '.mat'];
save(mat_3d_name,'ll_cube') ;

[mv,id] = max(ll_cube(:));
[i,j,d] = ind2sub(size(ll_cube),id)

surf(ll_cube(:,:,d))
% surf(squeeze(ll_cube(i,:,:)))
% surf(squeeze(ll_cube(:,j,:)))

parameters = [b_k(i) b_theta(j) b_wd(d)]

% k = find(ll_cube>=-1236400);
% [i,j,d] = ind2sub(size(ll_cube),k);
% ll_cube(i,j,d)=-2000000;
% [mv,id] = max(ll_cube(:));
% [i,j,d] = ind2sub(size(ll_cube),id)
% ll_cube(i,j,d)

% beta_0 = [0 1 2 0]
% beta_0 = [-100 100 -10]
% beta_0 = [-5.6814 2.6098 -0.0040]
% [b, hessian, grad, standard_error, covariance_matrix, t_stat, exit_flag, output] = band_runbi_ll_p(X, y, dummy_X, beta_0);
% [b, hessian, grad, standard_error, covariance_matrix, t_stat, exit_flag, output] = band_runbi_ll_p(X, y, beta_0);
% 
% save('b.mat','b') ;
% save('standard_error.mat','standard_error') ;
% save('t_stat.mat','t_stat') ;
% save('exit_flag.mat','exit_flag') ;
% 
% display(b)
% % display(standard_error)
% display(t_stat)
% display(grad)
% display(output)
% 
% m = grad'*(-inv(hessian))*grad;          % convergence criterion of Train

% diary off 

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
   
    