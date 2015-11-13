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
LR_test_result = 1-chi2cdf(likelihood_ratio_test,28)

% Or compare with the baseline model with 5 parameters including the constant
% ll_r = -286972                 % baseline mdoel (k = 5)
% ll_r = -286529                 
ll_r = -286471                    % model without three-way interaction terms (k = 25)
% ll_r = -286465
% ll_ur = -286459                  % full model
% ll_ur = -286461                     % model with 1st 3-way interation term (k = 26)
% ll_ur = -286467                     % model with 2nd 3-way interation term(k = 26)
% ll_ur = -286470                     % model with 3rd 3-way interation term(k = 26)
ll_ur = -286471                     % model with 4th 3-way interation term(k = 26)
% w: Note that the ll_ur results for 2-4th 3-way interaction terms were
% estimated using the starting point = 0 for the variable. This may influence the result
likelihood_ratio_test = -2*(ll_r-ll_ur)
chi2inv(0.99,4)
LR_test_result = 1-chi2cdf(likelihood_ratio_test,4)

% n = size(IVs,2)
n = 6885477
% k_ur = 29
k_ur = 26
k_r = 25
% k_r = 5
% ll = ll_ur
% ll = ll_r
AIC_r = -2*ll_r+2*k_r
BIC_r = -2*ll_r+ log(n)*k_r
AIC_ur = -2*ll_ur+2*k_ur
BIC_ur = -2*ll_ur+ log(n)*k_ur
% format long

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

% clear all
% load('innov_contin2.mat');
% load('explor_contin2.mat');

mean(innov_contin)
mean(explor_contin)
quantile(innov_contin(:,1),[0 0.05 0.1 0.25 0.3 0.5 0.6 0.75 0.9 0.95 1])
quantile(innov_contin(:,3),[0 0.05 0.1 0.25 0.3 0.5 0.6 0.75 0.9 0.95 1])
quantile(explor_contin(:,1),[0 0.05 0.1 0.25 0.3 0.5 0.6 0.75 0.9 0.95 1])
quantile(explor_contin(:,3),[0 0.05 0.1 0.25 0.3 0.5 0.6 0.75 0.9 0.95 1])
% quantile(innov_contin(:,1),[0.14 0.5 0.86])
% quantile(innov_contin(:,3),[0.14 0.5 0.86])
% quantile(explor_contin(:,1),[0.14 0.5 0.86])
% quantile(explor_contin(:,3),[0.14 0.5 0.86])

hist(innov_contin(:,1),30)
hist(innov_contin(:,3),30)
hist(explor_contin(:,1),30)
hist(explor_contin(:,3),30)
line([1.9 1.9],[0 16*10^5], 'Color', 'r','Linestyle',':')
line([-0.84 -0.84],[0 16*10^5], 'Color', 'r','Linestyle','--')

% WD_plug = mean_wd;
% WD_plug = 0;
base_prob_plug = mean_base_prob      %+1;
% innov_m_plug = mean_innov_m;
% innov_f_plug = mean_innov_f;
% explor_m_plug = mean_explor_m;
% explor_f_plug = mean_explor_f;
% innov_intv_m = [-0.8901   -0.2244    1.7034];     % 14% 50% 86% quantile
% innov_intv_f = [-0.9770   -0.1805    1.0119];     % 14% 50% 86% quantile
% explor_intv_m = [-0.8495   -0.4930    1.9075];    % 14% 50% 86% quantile
% explor_intv_f = [-0.7377   -0.2972    0.6773];    % 14% 50% 86% quantile
innov_intv_m = [-1   0    1];     % 14% 50% 86% quantile
innov_intv_f = [-1   0    1];     % 14% 50% 86% quantile
explor_intv_m = [-1   0    1];    % 14% 50% 86% quantile
explor_intv_f = [-1   0    1];    % 14% 50% 86% quantile

% WD_plug = -10:1:40;
WD_plug = 2;

prob = zeros(length(WD_plug),length(innov_intv_m));
for q = 1:length(innov_intv_m)
    innov_m_plug = innov_intv_m(2);
    innov_f_plug = innov_intv_f(2);
%     explor_m_plug = explor_intv_m(q);
%     explor_f_plug = explor_intv_f(q);
    explor_m_plug = explor_intv_m(2);            % choose median score for other variables
    explor_f_plug = explor_intv_f(q);
    
    innov_plug = [innov_m_plug innov_m_plug^2 innov_f_plug innov_f_plug^2 innov_m_plug*innov_f_plug (innov_m_plug*innov_f_plug)^2];
    explor_plug = [explor_m_plug explor_m_plug^2 explor_f_plug explor_f_plug^2 explor_m_plug*explor_f_plug (explor_m_plug*explor_f_plug)^2];
    
    for k = 1:length(WD_plug)
        week_IV = 100*gampdf(WD_plug(k),exp(b(2)),exp(b(3)));
        week_IV(WD_plug(k)<1)=0;
        
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
        prob(k,q)=1./(1+exp_util);
        % this is still the probability of choosing the product
        % pmat = [prob 1-prob];
        % pmat = pmat.*choice_dv;
        % [r c p] = find(pmat);                                             % I*1
        % ll_ur = sum(log(p))                 % ll of unrestricted model (the one that was estimated)
    end
end
plot(WD_plug, prob(:,1),'r')
hold on 
plot(WD_plug, prob(:,2),'g')
plot(WD_plug, prob(:,3),'b')
legend('innov_m = L','innov_m = M','innov_m = H')
title('predicted likelihood plot (friend innovativeness is low)')
xlabel('week diff')
ylabel('probability')
set(gca, 'YLim',[0 0.014])                    % change y axis
hold off

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

%%
age = 0:1:70;
male = 1;
w = -3+0.1*age+0.5*male;
y = exp(w)./(1+exp(w));
plot(age,y, 'r')
hold on
male = 0;
w = -3+0.1*age+0.5*male;
y = exp(w)./(1+exp(w));
plot(age,y, 'b')
hold off

w = -5:0.1:5;
y = exp(w)./(1+exp(w));
plot(y,'r')
hold on
w = w+1;
y = exp(w)./(1+exp(w));
plot(y,'b')
hold off

%% corr time appear in sample vs explor_score
clear all

load('matpstd2.mat');
row_ind = 1
row_fid = [];
row_interval = [];
new_row = matp(1,2);
for i = 2:size(matp,1)
    old_row = new_row;
    new_row = matp(i,2);
    if new_row == old_row
        row_ind = row_ind + 1;
    else
        row_fid = [row_fid old_row];
        row_interval = [row_interval row_ind];
        row_ind = 1;
    end
end
row_fid = [row_fid new_row];
row_interval = [row_interval row_ind];

row_interval_sum = cumsum(row_interval);

save('row_fid.mat','row_fid');
save('row_interval.mat','row_interval');

%
row_ind = 1
row_mid = [];
row_num = [];
new_row = matp(1,1);
for i = 2:size(matp,1)
    old_row = new_row;
    new_row = matp(i,1);
    if new_row == old_row
        row_ind = row_ind + 1;
    else
        row_mid = [row_mid old_row];
        row_num = [row_num row_ind];
        row_ind = 1;
    end
end
row_mid = [row_mid new_row];
row_num = [row_num row_ind];

sum(row_num)

row_num_sum = cumsum(row_num);

save('row_num.mat','row_num');
save('row_mid.mat','row_mid');
% csvwrite('row_num.csv',row_num);
% csvwrite('row_mid.csv',row_mid);

load('innov_contin_std2.mat');
load('explor_contin_std2.mat');

innov_m_score = innov_contin(row_num_sum,1);
innov_f_score = innov_contin(row_interval_sum,3);
explor_m_score = explor_contin(row_num_sum,1);
explor_f_score = explor_contin(row_interval_sum,3);
innov_m_score2 = innov_contin(row_interval_sum,1).^2;
explor_m_score2 = explor_contin(row_interval_sum,1).^2;

corr(innov_m_score,row_num')
corr(innov_f_score,row_interval')
corr(explor_m_score,row_num')
corr(explor_f_score,row_interval')
corr(innov_m_score2,row_interval')
corr(explor_m_score2,row_interval')

corr(innov_m_score,explor_m_score)            % 0.03
corr(innov_f_score,explor_f_score)            % 0.10
% compare with
innov = csvread('EAi3.csv');
explor = csvread('explorer3.csv');
corr(innov(:,2), explor)                      % -0.04

% quantile(innov_f_score,[0 0.05 0.25 0.5 0.75 0.95 1])
% quantile(explor_f_score,[0 0.05 0.25 0.5 0.75 0.95 1])
% 
% quantile(innov_m_score,[0 0.05 0.25 0.5 0.75 0.95 1])
% quantile(explor_m_score,[0 0.05 0.25 0.5 0.75 0.95 1])
% 
% corr(matp(:,5),innov_contin(:,3))
% corr(matp(:,5),explor_contin(:,3))
% corr(matp(:,5),innov_contin(:,1))
% corr(matp(:,5),explor_contin(:,1))

corr(innov_contin(row_interval_sum,1),explor_contin(row_interval_sum,1))         
corr(innov_contin(row_interval_sum,2),explor_contin(row_interval_sum,2))
% high corr bewteen member innov and explor scores may explain why the
% main effect cofficients are insig
% corr(innov_contin(:,1),explor_contin(:,1))             
% corr(innov_contin(:,2),explor_contin(:,2))
% corr(innov_contin(:,3),explor_contin(:,3))
% corr(innov_contin(:,4),explor_contin(:,4))