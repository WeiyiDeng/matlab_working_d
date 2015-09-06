clc
clear all

load('matp.mat');
% matlabpool local 2

X = matp(:,[6 7 8]);
y = matp(:,5);

IVs = X;
choice_dv = [y 1-y];

plus = zeros(length(choice_dv),1);
plus(IVs(:,2)==0)=0.01;
IVs(:,2) = IVs(:,2)+plus;
clearvars plus
% [row col val] = find(IVs(:,2)==0);
% IVs(row,2) = 0.01;

tic
b = [-6.1474 0.6 4 3.1975 2];

const = b(1);
bs = b(4:end)';

% interval = length(IVs(:,2))/num_workers;
% ind = interval;

week_IV = IVs(:,3).*gampdf(IVs(:,2),b(2),b(3));
% week_IV = zeros(IVs(:,2),1);
% parfor i = 1:length(IVs(:,2))
%     week_IV(i) = IVs(i,3)*gampdf(IVs(i,2),b(2),b(3));
% end
    
FV = [IVs(:,1) week_IV]*bs;

exp_util = exp(-(const+FV));         % this is now the utility of the external good
prob=1./(1+exp_util);                % this is still the probability of choosing the product
pmat = [prob 1-prob];
pmat = pmat.*choice_dv;
[r c p] = find(pmat);                                             % I*1
ll = sum(log(p));
toc

% sum(isnan(week_IV))
% ind = find(isnan(week_IV));
% values = IVs(ind,2);
% 
% gampdf(IVs(ind(1),2),0.6,4)
