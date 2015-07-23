clc
clear all

% load('mat1.mat');
% load('mat2.mat');

band_mat1 = csvread('mat1.csv');
band_mat2 = csvread('mat2.csv');

member_inds1 = csvread('member_ind_mat1.csv');
member_inds2 = csvread('member_ind_mat2.csv');

% 153-1 members left in mat1 and mat2
% member 2978 in mat1 has all 0 adoptions DV (with friend 3789 and 1026)
% both of them did not adopt any song during the member-friend overlap of 
% the second half of their active periods
% some members in mat2 have less than 30 nobs, removed from mat2
% 144 members left in mat2mod.csv 

%%
% I = 152
I = 3
K = 2
b_store = zeros(I,K+1);
se_store = zeros(I,K+1);
cov_store = zeros(K+1,K+1,I);
tstat_store = zeros(I,K+1);
hess = zeros(K+1,K+1,I);

i_id = zeros(I,1);
% beta_0 = zeros(1,K+1);
beta_0 = [0 0 0];
ind = 0;
for i = 1:I
    i_id(i) = band_mat1(ind+1,1);
    nrows = member_inds1(member_inds1(:,1)==i_id(i),2);
    imat = band_mat1(ind+1:ind+nrows,:);         % ind tbc
%     csvwrite('imat1.csv',imat);
    ind = ind+nrows;
    X = imat(:,[5 7]);
    y = imat(:,4);
    [b, hessian, standard_error, covariance_matrix, t_stat] = band_runbi_ll(X, y, beta_0);
    b_store(i,:) = b;
    se_store(i,:) = standard_error';
    cov_store(:,:,i) = covariance_matrix;
    tstat_store(i,:) = t_stat;
    hess(:,:,i) = hessian;
end

% [x,fval,exitflag,output,grad,hessian] = logregr(imat(:,5), imat(:,7), imat(:,4), 0, 0, 0)
% [b,dev,stats] = glmfit([imat(:,5) imat(:,7)], imat(:,4),'binomial')