clc
% clear
clear all

% diary('wwdiary.txt')

% date_=clock;
% resultsfilename=['Results/r_Results_' num2str(date_(1)) '_' num2str(date_(2)) '_' num2str(date_(3)) '-' num2str(date_(4))  '_' num2str(date_(5)) '.txt'];
% diary(resultsfilename);

%% goodness of fit
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

McFadden_R_square = 1 - ll_ur/ll_r                        % likelihood ratio index
% pseudo_R_square = 1 - 1/(1+2*(ll_ur-ll_r)/length(y))

% format long g

likelihood_ratio_test = -2*(ll_r-ll_ur)                   % >> chi2inv(0.95,28)

% Or compare with the baseline model with 5 parameters including the constant
% ll_r = -286972                 % baseline mdoel (k = 5)
% ll_r = -286529                 
ll_r = -286471                   % model without three-way interaction terms (k = 25)
ll_ur = -286459                  % full model
likelihood_ratio_test = -2*(ll_r-ll_ur)
chi2inv(0.99,4)

% n = size(IVs,2)
n = 6885477
k = 29
% k = 25
ll = ll_ur
% ll = ll_r
AIC = -2*ll+2*k
BIC = -2*ll+ log(n)*k

% Percent of correctly predicted values not reported because it is against
% the meaning of specifying choice probabilities (researcher cannot predict
% the choice decision but only the probabilities of each choice) see Train Po69 

% A = prob > 0.5;
% sum(A)

% marginal effect of baseline probability to adopt band j at time t
exp_util2 = exp(const+FV);         
prob2 = exp_util2./(1+exp_util2);               

marginal_effect_based_prob = mean(exp_util2./(1+exp_util2).^2)*b(4)
% One std change in baseline probability the               probability of
% adopting band j at time t increases 0.0023 (from probability such as 0.0061)
%% coefficients
clear all

load('matpstd2.mat');
load('innov_contin_std2.mat');
load('explor_contin_std2.mat');

% sum(innov_contin(:,1)^2~=innov_contin(:,2))           % testing

mean_base_prob = mean(matp(:,6));                       % sd = 1
mean_wd = mean(matp(:,7))                              
std_wd = std(matp(:,7));
mean_innov_m = mean(innov_contin(:,1));                 % sd = 1
mean_innov_f = mean(innov_contin(:,3));                 % sd = 1
mean_explor_m = mean(explor_contin(:,1));               % sd = 1
mean_explor_f = mean(explor_contin(:,3));               % sd = 1

b = [-5.0951    0.19983    2.0459    0.3222    0.0627];
b = [b     -0.0102    0.0667   -0.1275    0.0555   -0.0005   -0.0652    -0.0225    0.0052];  
b = [b      0.0181   -0.0185    0.0131   -0.0125   -0.0088    0.0092   -0.0145    0.0017];
b = [b     -0.0160   -0.0046    0.0027    0.0005    0.0117    0.0008    0.0052   -0.0014]

const = b(1);

b_basic = b(4:5)';
b_innov = b([6:9 22:23 14:17 26:27])';
b_explor = b([10:13 24:25 18:21 28:29])';

quantile(explor_contin(:,3),[0,0.16,0.84,0.975,1])
% A = quantile(explor_contin(:,3),0.25)

% WD_plug = mean_wd;
WD_plug = 2;
base_prob_plug = mean_base_prob      %+1;
innov_m_plug = mean_innov_m;
innov_f_plug = mean_innov_f;
explor_m_plug = mean_explor_m;
explor_f_plug = mean_explor_f;

innov_plug = [innov_m_plug innov_m_plug^2 innov_f_plug innov_f_plug^2 innov_m_plug*innov_f_plug (innov_m_plug*innov_f_plug)^2];
explor_plug = [explor_m_plug explor_m_plug^2 explor_f_plug explor_f_plug^2 explor_m_plug*explor_f_plug (explor_m_plug*explor_f_plug)^2];
week_IV = 100*gampdf(WD_plug,exp(b(2)),exp(b(3)));          
week_IV(WD_plug<1)=0;           

innov_WD_multip = zeros(size(innov_plug));
for i = 1:size(innov_plug,2);
    innov_WD_multip(:,i) = innov_plug(:,i).*week_IV;
end

explor_WD_multip = zeros(size(explor_plug));
for j = 1:size(explor_plug,2);
    explor_WD_multip(:,j) = explor_plug(:,j).*week_IV;
end

FV = [base_prob_plug week_IV]*b_basic + innov_plug*b_innov(1:6)+innov_WD_multip*b_innov(7:end) + explor_plug*b_explor(1:6)+explor_WD_multip*b_explor(7:end);

% % another way of calculating marginal effect of baseline probability
% exp_util = exp(const+FV);         
% prob = exp_util/(1+exp_util)               
% 
% marginal_effect_based_prob = exp_util/(1+exp_util)^2*b(4)

exp_util = exp(-(const+FV));         % this is now the utility of the external good
prob=1./(1+exp_util)                % this is still the probability of choosing the product
% pmat = [prob 1-prob]; 
% pmat = pmat.*choice_dv;
% [r c p] = find(pmat);                                             % I*1
% ll_ur = sum(log(p))                 % ll of unrestricted model (the one that was estimated)


%% analytical derivative
syms X Z
w = 3*X+4*Z+5*X^2*Z^2
dP_dX = diff(w,X)
dP_dXdZ = diff(dP_dX,Z)
deriv_f = matlabFunction(dP_dXdZ)
deriv_f(2,3)

b = [-5.0951    0.19983    2.0459    0.3222    0.0627];
b = [b     -0.0102    0.0667   -0.1275    0.0555   -0.0005   -0.0652    -0.0225    0.0052];  
b = [b      0.0181   -0.0185    0.0131   -0.0125   -0.0088    0.0092   -0.0145    0.0017];
b = [b     -0.0160   -0.0046    0.0027    0.0005    0.0117    0.0008    0.0052   -0.0014]

const = b(1);

% b_basic = b(4:5)';
% b_innov = b([6:9 22:23 14:17 26:27])';
% b_explor = b([10:13 24:25 18:21 28:29])';

base_prob = 0.005;

syms Xim Xif WD

% WD = 100*gampdf(44,exp(b(2)),exp(b(3)));          
% WD(44<1)=0;

Xem = 0.04
Xef = 0.3

% innov_X = [Xim Xim^2 Xif Xif^2 Xim*Xif (Xim*Xif)^2];
% explor_X = [Xem Xem^2 Xef Xef^2 Xem*Xef (Xem*Xef)^2];
% 
% 
% innov_X = [Xim Xim^2 Xif Xif^2 Xim*Xif (Xim*Xif)^2];
% explor_X = [Xem Xem^2 Xef Xef^2 Xem*Xef (Xem*Xef)^2];
% 
% innov_WD_multip = innov_X.*WD;
% 
% explor_WD_multip = explor_X.*WD;

FV = b(4)*base_prob + b(5)*WD + b(6)*Xim+ b(7)*Xim^2+ b(8)*Xif+ b(9)*Xif^2+ b(22)*Xim*Xif+ b(23)*(Xim*Xif)^2+ b(14)*Xim*WD+ b(15)*Xim^2*WD+ b(16)*Xif*WD+ b(17)*Xif^2*WD+ b(26)*Xim*Xif*WD+ b(27)*(Xim*Xif)^2*WD+ b(10)*Xem+ b(11)*Xem^2+ b(12)*Xef+ b(13)*Xef^2+ b(24)*Xem*Xef+ b(25)*(Xem*Xef)^2+ b(18)*Xem*WD+ b(19)*Xem^2*WD+ b(20)*Xef*WD+ b(21)*Xef^2*WD+ b(28)*Xem*Xef*WD+ b(29)*(Xem*Xef)^2*WD

% FV = [base_prob WD]*b_basic + innov_X*b_innov(1:6)+innov_WD_multip*b_innov(7:end) + explor_X*b_explor(1:6)+explor_WD_multip*b_explor(7:end)

w = const + FV
prob = exp(w)/(1+exp(w))
dP_dXim = diff(prob,Xim)
dP_dXimdXif = diff(dP_dXim,Xif)
deriv_f = matlabFunction(dP_dXimdXif)
deriv_f(2,0,0)




