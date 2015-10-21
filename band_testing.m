clc
% clear
clear all

% diary('wwdiary.txt')

% date_=clock;
% resultsfilename=['Results/r_Results_' num2str(date_(1)) '_' num2str(date_(2)) '_' num2str(date_(3)) '-' num2str(date_(4))  '_' num2str(date_(5)) '.txt'];
% diary(resultsfilename);

% load('matp_b.mat');
load('matpstd2.mat');

% X = matp(:,[6 7 8]);
y = matp(:,5);

% IVs = X;
choice_dv = [y 1-y];

b = -4.9382;

const = repmat(b,length(y),1);

exp_util = exp(-const);         % this is now the utility of the external good
prob=1./(1+exp_util);                % this is still the probability of choosing the product
pmat = [prob 1-prob];
pmat = pmat.*choice_dv;
[r c p] = find(pmat);                                             % I*1
ll_r = sum(log(p))              % restricted ll of model with only the intercept (all other parameters set to 0)

% ur
X = matp(:,[6 7]);
y = matp(:,5);

IVs = X;
choice_dv = [y 1-y];

clearvars clearvars -EXCEPT ll_r matp

load('innov_contin_std2.mat');
load('explor_contin_std2.mat');

innov_X = innov_contin;
explor_X = explor_contin;

clearvars innov_contin explor_contin

b = [-5.0951    0.19983    2.0459    0.3222    0.0627];
b = [b     -0.0102    0.0667   -0.1275    0.0555   -0.0005   -0.0652    -0.0225    0.0052];  
b = [b      0.0181   -0.0185    0.0131   -0.0125   -0.0088    0.0092   -0.0145    0.0017];
b = [b     -0.0160   -0.0046    0.0027    0.0005    0.0117    0.0008    0.0052   -0.0014]

const = b(1);

b_basic = b(4:5)';
b_innov = b([6:9 22:23 14:17 26:27])';
b_explor = b([10:13 24:25 18:21 28:29])';

week_IV = 100*gampdf(IVs(:,2),exp(b(2)),exp(b(3)));          
week_IV(IVs(:,2)<1)=0;

innov_WD_multip = zeros(size(innov_X));
for i = 1:size(innov_X,2);
    innov_WD_multip(:,i) = innov_X(:,i).*week_IV;
end

explor_WD_multip = zeros(size(explor_X));
for j = 1:size(explor_X,2);
    explor_WD_multip(:,j) = explor_X(:,j).*week_IV;
end

FV = [IVs(:,1) week_IV]*b_basic + innov_X*b_innov(1:6)+innov_WD_multip*b_innov(7:end) + explor_X*b_explor(1:6)+explor_WD_multip*b_explor(7:end);

exp_util = exp(-(const+FV));         % this is now the utility of the external good
prob=1./(1+exp_util);                % this is still the probability of choosing the product
pmat = [prob 1-prob]; 
pmat = pmat.*choice_dv;
[r c p] = find(pmat);                                             % I*1
ll_ur = sum(log(p))                 % ll of unrestricted model (the one that was estimated)

McFadden_R_square = 1 - ll_ur/ll_r