clc
clear all

% load('mat1b.mat');
% load('mat2b.mat');
% load('member_adopt_b.mat');        % w: error message cannot read file??

% mat1 = csvread('mat1.csv');

mat1 = csvread('mat1.csv');
mat2 = csvread('mat2.csv');      % 13137715	13145733

% mat1(13137715:13145733,:)=[];           % removes obs of member 2978
% member_adopt(13137715:13145733,:)=[];

member_inds1 = csvread('member_ind_mat1.csv');
member_inds2 = csvread('member_ind_mat2.csv');

% 153-1 members left in mat1 and mat2
% member 2978 in mat1 has all 0 adoptions DV (with friend 3789 and 1026)
% both of them did not adopt any song during the member-friend overlap of 
% the second half of their active periods
% some members in mat2 have less than 30 nobs, removed from mat2
% 144 members left in mat2mod.csv 

% mat1_dummies = [EAm_a.*MAf_a EAm_a.*LAf_a MAm_a.*EAf_a MAm_a.*MAf_a MAm_a.*LAf_a LAm_a.*EAf_a LAm_a.*MAf_a LAm_a.*LAf_a];
% mat1_dummies = [mat1(:,8).*mat1(:,12) mat1(:,8).*mat1(:,13) mat1(:,9).*mat1(:,11) mat1(:,9).*mat1(:,12) mat1(:,9).*mat1(:,13) mat1(:,10).*mat1(:,11) mat1(:,10).*mat1(:,12) mat1(:,10).*mat1(:,13)];
% band_mat1 = [mat1(:,1:7) mat1_dummies];
% clearvars mat1

%%
% I = 152
NUSER = 8320;                % max number of users
NBAND = 6046;                % max number of bands
% I = 152
I = 3
% K = 4                        % number of parameters not including constant
% b_store = zeros(I,K+1);
% se_store = zeros(I,K+1);
% cov_store = zeros(K+1,K+1,I);
% tstat_store = zeros(I,K+1);
% hess = zeros(K+1,K+1,I);
b_store = cell(I,1);
se_store = cell(I,1);
cov_store = cell(I,1);
tstat_store = cell(I,1);
hess = cell(I,1);
residuals = zeros(size(mat1,1),1);
resid_mat = cell(NUSER,1);
EA_location = cell(I,1);
f_types = cell(I,1);
m_type = cell(I,1);
% cell_spacealloc = max(member_inds1(:,2));
% cell_spacealloc = ceil(quantile(member_inds1(:,2),0.8));

i_id = zeros(I,1);
% beta_0 = zeros(1,K+1);
est1_exitflag = zeros(I,1);
% beta_0 = [0 0 0];
ind = 0;
for i = 1:I
    i_id(i) = mat1(ind+1,1);
    nrows = member_inds1(member_inds1(:,1)==i_id(i),2);
    imat = mat1(ind+1:ind+nrows,:);         % ind tbc
    m_type{i} = imat(1,8:10);
    f_EA = any(imat(:,11)) == 1;
    f_MA = any(imat(:,12)) == 1;
    f_LA = any(imat(:,13)) == 1;
    f_types{i} = [f_EA f_MA f_LA];
    if f_EA == 1 && f_MA == 1 && f_LA == 1
        EA_col = [11:12];
    elseif f_EA == 1 && f_MA == 1
        EA_col = 11;
    elseif f_MA == 1 && f_LA == 1
        EA_col = 12;
    elseif f_EA == 1 && f_LA == 1
        EA_col = 11;
    else
        EA_col = [];
    end
    EA_location{i} = EA_col;

%     if any(imat(:,11)) == 1 && any(imat(:,12)) == 1 && any(imat(:,13)) == 1
%         EA_col = [11:12];
%     elseif any(imat(:,11)) == 1 && any(imat(:,12)) == 1
%         EA_col = 11;
%     elseif any(imat(:,12)) == 1 && any(imat(:,13)) == 1
%         EA_col = 12;
%     elseif any(imat(:,11)) == 1 && any(imat(:,13)) == 1
%         EA_col = 11;
%     else
%         EA_col = [];
%     end
%     EA_location{i} = EA_col;
    %     csvwrite('imat1.csv',imat);
%     X = imat(:,[5 7 11:12]);
%     X = [member_adopt(ind+1:ind+nrows,:) imat(:,[5 7 EA_col])];
    X = imat(:,[5 7 EA_col]);
    y = imat(:,4);
    beta_0 = zeros(1,size(X,2)+1);
    [b, hessian, standard_error, covariance_matrix, t_stat, exit_flag] = band_runbi_ll(X, y, beta_0);
%     b_store(i,:) = b;
%     se_store(i,:) = standard_error';
%     cov_store(:,:,i) = covariance_matrix;
%     tstat_store(i,:) = t_stat;
%     hess(:,:,i) = hessian;
%     est1_exitflag(i) = exit_flag;
    b_store{i} = b;
    se_store{i} = standard_error';
    cov_store{i} = covariance_matrix;
    tstat_store{i} = t_stat;
    hess{i} = hessian;
    est1_exitflag(i) = exit_flag;
    util_hat = b(1)+ X*b(2:end)';
    %     yhat = exp(util_hat)./(1+exp(util_hat));
    yhat = 1./(1+exp(-util_hat));
    resid = y-yhat;
    residuals(ind+1:ind+nrows) = resid;
%     resid_mat{i_id(i)} = sparse(imat(:,2),imat(:,3),resid,NUSER,NBAND,cell_spacealloc);
    resid_mat{i_id(i)} = sparse(imat(:,2),imat(:,3),resid,NUSER,NBAND);
    ind = ind+nrows;
    if i >= 10
        break
    end
end

% [x,fval,exitflag,output,grad,hessian] = logregr(imat(:,5), imat(:,7), imat(:,4), 0, 0, 0)
% [b,dev,stats] = glmfit([imat(:,5) imat(:,7)], imat(:,4),'binomial')

save('resid_mat.mat','resid_mat') ;