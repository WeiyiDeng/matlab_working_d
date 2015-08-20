clc
clear all

load('matp.mat');
% % mat1 = csvread('imat4474.csv');
% 
% X = matp(:,[6 7]);
% y = matp(:,5);
% % X = mat1(:,[5 7]);
% % y = mat1(:,4);
% beta_0 = zeros(1,size(X,2)+1);
% [b, hessian, standard_error, covariance_matrix, t_stat, exit_flag] = band_runbi_ll(X, y, beta_0);
% 
% save('b.mat','b') ;
% save('standard_error.mat','standard_error') ;
% save('t_stat.mat','t_stat') ;
% save('exit_flag.mat','exit_flag') ;

%%
load('row_mid.mat');
load('row_num.mat');

I = 158

b_store_p = cell(I,1);
se_store_p = cell(I,1);
cov_store_p = cell(I,1);
tstat_store_p = cell(I,1);
% hess = cell(I,1);
i_id = zeros(I,1);
est_exitflag_p = zeros(I,1);

ind = 0;
for i = 1:I
    i_id(i) = row_mid(i);
    nrows = row_num(i);
    imat = matp(ind+1:ind+nrows,:);         % ind tbc
    
    X = imat(:,[6 7]);
    y = imat(:,5);
    beta_0 = zeros(1,size(X,2)+1);
    [b, hessian, standard_error, covariance_matrix, t_stat, exit_flag] = band_runbi_ll(X, y, beta_0);

    b_store_p{i} = b;
    se_store_p{i} = standard_error';
    cov_store_p{i} = covariance_matrix;
    tstat_store_p{i} = t_stat;
%    hess{i} = hessian;
    est_exitflag_p(i) = exit_flag;
        
    ind = ind+nrows; 
%     if i >= 10
%         break
%     end    
end
   
save('b_store_p.mat','b_store_p') ;
save('se_store_p.mat','se_store_p') ;
save('tstat_store_p.mat','tstat_store_p') ;
save('est_exitflag_p.mat','est_exitflag_p') ;
   
    