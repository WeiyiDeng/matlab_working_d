%%
% h = fopen('memmapex.raw','wb');
% for ii = 1:500
% fwrite(h,randn(60*256*256,1),'double');
% end
% fclose(h);
% 
% idx = randperm(30000);
% idx = sort(idx(1:3000));
%  
% rmean = zeros(256,256);
% m = memmapfile('memmapex.raw','Format',{'double',[256,256],'im'});
%  
% for ii = 1:length(idx)
% rmean = rmean + m.Data(idx(ii)).im;
% if mod(ii,100) == 0
% ii
% end
% end
%  
% rmean = rmean/length(idx);

%%
% % load('innov_contin_std.mat')
% % fileID = fopen('innov_contin_stdat.dat','w');
% % fwrite(fileID, innov_contin,'double');
% % fclose(fileID);
% % clear('fileID')
% % 
% % load('explor_contin_std.mat')
% % fileID = fopen('explor_contin_stdat.dat','w');
% % fwrite(fileID, explor_contin,'double');
% % fclose(fileID);
% % clear('fileID')
% 
% clear all
% load('matpstd.mat');
% X = matp(:,[6 7]);
% y = matp(:,5);
% 
% clearvars matp
% 
% tic
% IVs = X;
% choice_dv = [y 1-y];
% 
% beta_0 = [-5.6271    0.8375   27.4694    0.3870    0.1020];
% beta_0 = [beta_0         -0.0509    0.0281   -0.2486    0.2183    0.1682   -0.1598    -0.0583    0.0421];  
% beta_0 = [beta_0          0.0102    0.0140    0.0414   -0.0358    0.0538   -0.0406    0.0038   -0.0118];
% beta_0 = [beta_0         -0.0112   -0.0003   -0.0245    0.0172    0.0076    0.0081   -0.0026    0.0057]; 
% b = beta_0; 
% 
% const = b(1);
% 
% b_basic = b(4:5)';
% b_innov = b([6:9 22:23 14:17 26:27])';
% b_explor = b([10:13 24:25 18:21 28:29])';
% 
% week_IV = 100*gampdf(IVs(:,2),b(2),b(3));          
% week_IV(IVs(:,2)<1)=0;
% 
% FV_basic = [IVs(:,1) week_IV]*b_basic;
% 
% % load('innov_contin_std.mat');
% mm = memmapfile('innov_contin_stdat.dat', 'Format', {'double', [47085403 6], 'col'},'Repeat', 1);
% % m = mm.Data(1).col;          % first col
% % info = whos('mm');
% m = mm.Data;
% 
% innov_WD_multip = zeros(size(m.col));
% for j = 1:6;
%     innov_WD_multip(:,j) = m.col(:,j).*week_IV;
% end
% % innov_WD_multip = [];
% 
% FV_innov = m.col*b_innov(1:6)+innov_WD_multip*b_innov(7:end);
% 
% clearvars innov_contin innov_WD_multip m mm
% 
% % load('explor_contin_std.mat');
% % explor_X = explor_contin(:,1:4);
% mm = memmapfile('explor_contin_stdat.dat', 'Format', {'double', [47085403 6], 'col'},'Repeat', 1);
% % m = mm.Data(1).col;          % first col
% % info = whos('mm');
% m = mm.Data;
% 
% explor_WD_multip = zeros(size(m.col));
% for j = 1:6;
%     explor_WD_multip(:,j) = m.col(:,j).*week_IV;
% end
% % explor_WD_multip = [];
% 
% FV_explor = m.col*b_explor(1:6)+explor_WD_multip*b_explor(7:end);
% 
% clearvars explor_contin explor_WD_multip m mm
% 
% FV = FV_basic + FV_innov + FV_explor;
% 
% % exp_util = exp(const+FV);          % utility of choosing the product
% % prob=exp_util./(1+exp_util);
% exp_util = exp(-(const+FV));         % this is now the utility of the external good
% prob=1./(1+exp_util);                % this is still the probability of choosing the product
% pmat = [prob 1-prob]; 
% pmat = pmat.*choice_dv;
% [r c p] = find(pmat);                                             % I*1
% LL = -sum(log(p));     
% 
% toc

%%
clear all
load('matpstd.mat');
X = matp(:,[6 7]);
y = matp(:,5);

nrows = size(matp,1);

clearvars matp

% load('innov_contin_std.mat');
% 
% fileID = fopen('mytry1.dat','w');
% fwrite(fileID, matp,'double');
% fclose(fileID);
% 
% mm = memmapfile('mytry1.dat', 'Format', {'double', [nrows 1], 'mj'},'Repeat', 6);
% A = mm.Data(1).mj;          % first col
% dtmean = [];
% for i = 1:8
%     dtmean = [dtmean mean(mm.Data(i).mj)];
% end

tic
IVs = X;
choice_dv = [y 1-y];

beta_0 = [-5.6271    0.8375   27.4694    0.3870    0.1020];
beta_0 = [beta_0         -0.0509    0.0281   -0.2486    0.2183    0.1682   -0.1598    -0.0583    0.0421];  
beta_0 = [beta_0          0.0102    0.0140    0.0414   -0.0358    0.0538   -0.0406    0.0038   -0.0118];
beta_0 = [beta_0         -0.0112   -0.0003   -0.0245    0.0172    0.0076    0.0081   -0.0026    0.0057]; 
b = beta_0 

const = b(1);

b_basic = b(4:5)';
b_innov = b([6:9 22:23 14:17 26:27])';
b_explor = b([10:13 24:25 18:21 28:29])';

week_IV = 100*gampdf(IVs(:,2),b(2),b(3));          
week_IV(IVs(:,2)<1)=0;

FV_basic = [IVs(:,1) week_IV]*b_basic;

load('innov_contin_std.mat');

innov_WD_multip = zeros(size(innov_contin));
for i = 1:size(innov_contin,2);
    innov_WD_multip(:,i) = innov_contin(:,i).*week_IV;
end
% innov_WD_multip = [];

% FV_innov = [innov_contin innov_WD_multip]*b_innov;
FV_innov = innov_contin*b_innov(1:6)+innov_WD_multip*b_innov(7:end);

clearvars innov_contin innov_WD_multip

load('explor_contin_std.mat');
% explor_X = explor_contin(:,1:4);

explor_WD_multip = zeros(size(explor_contin));
for j = 1:size(explor_contin,2);
    explor_WD_multip(:,j) = explor_contin(:,j).*week_IV;
end
% explor_WD_multip = [];

% FV_explor = [explor_contin explor_WD_multip]*b_explor;
FV_explor = explor_contin*b_explor(1:6)+explor_WD_multip*b_explor(7:end);

clearvars explor_contin explor_WD_multip

FV = FV_basic + FV_innov + FV_explor;

% exp_util = exp(const+FV);          % utility of choosing the product
% prob=exp_util./(1+exp_util);
exp_util = exp(-(const+FV));         % this is now the utility of the external good
prob=1./(1+exp_util);                % this is still the probability of choosing the product
pmat = [prob 1-prob]; 
pmat = pmat.*choice_dv;
[r c p] = find(pmat);                                             % I*1
LL = -sum(log(p));     

toc
