% function [b, hessian, grad, standard_error, covariance_matrix, t_stat, exit_flag, output] = band_runbi_ll_p(X, y, dummy_X, beta_0)
function [b,fval,exit_flag,output,beta_hist] = band_runbi_ll_SA(X, y, beta_0, lb, ub)
global I J choice_dv IVs dummies se T b_best ll_best

% I = 10000
% J = 2
% T = 1
% % CST = ones(100,1).*T;
% 
% rng(0);
% 
% IV1 = randn(I*T,1);                             
% IV2 = randn(I*T,1);
% IVs = [IV1 IV2];
IVs = X;

choice_dv = [y 1-y];

plus = zeros(length(choice_dv),1);
plus(IVs(:,2)==0)=0.01;
IVs(:,2) = IVs(:,2)+plus;
clearvars plus

% dummies = dummy_X;
% b0 = 7;
% % b0 = repmat(b0,I*T,1);
% 
% b1 = 34;
% b2 = -21;
% fbs = [b1; b2];

% exp_util=exp(b0+ IVs*fbs);
% prob=1./(1+exp_util);
% draw_for_choice=rand(I*T,1);
% choice=prob<draw_for_choice;

% DV = choice;
% choice_dv = [DV 1-DV];

% beta0 = [0 0 0]
beta0 = beta_0;

% options = optimset('LargeScale','on','GradObj','on','Hessian','on','TolFun',1e-6, 'MaxIter',1e4, 'MaxFunEvals', 1e5)         % LargeScale off is quasi-Newton method in optimset
% options = optimset('LargeScale','off','GradObj','off','Hessian','off')
% options = optimset('Display','iter','LargeScale','off','GradObj','off','Hessian','off','TolFun',1e-6, 'TolX',1e-6, 'MaxIter',1e3, 'MaxFunEvals', 1e5, 'PlotFcns',@optimplotfirstorderopt);   %, 'FinDiffType', 'central');
% options = optimset('Display','iter','LargeScale','off','GradObj','off','Hessian','off','TolFun',1e-6, 'TolX',1e-6, 'MaxIter',0, 'MaxFunEvals', 1e5, 'PlotFcns',@optimplotfirstorderopt);   %, 'FinDiffType', 'central');

% saoptimset('simulannealbnd')         % w: check out the default options
b_best = beta0;
ll_best = 1.0000e+10;
save('b_best.mat','b_best')
save('ll_best.mat','ll_best')

% options = saoptimset('Display','iter','DisplayInterval',400,'OutputFcns',@(x,optimValues,state)beta_iter(x,optimValues,state));
options = saoptimset('Display','iter','DisplayInterval',400,'OutputFcns',@beta_iter);

% [b, fval,exitflag,output,grad,hessian] = fminunc(@band_bi_ll_p,beta0,options);
[b, fval,exitflag,output] = simulannealbnd(@band_bi_ll_p,beta0,lb,ub,options);
% [b, fval,exitflag,output] = simulannealbnd(@band_bi_ll_p,beta0,lb,ub);

% disp(['constant ' num2str(b(1)) '']);
% disp(['coefficients ' num2str(b(2:end)) '']);

% standard_error = sqrt(diag(inv(hessian)));
% % These s.e. are correct only when grad and hessian are NOT provided by me 
% % (turn off GradObj and Hessian to use these standard errors)
% % w: why is the case ??
% 
% covariance_matrix = inv(hessian);
% 
% t_stat = b./standard_error';
% disp(['t statistics ' num2str(t_stat) '']);

exit_flag = exitflag;

beta_hist = b_hist;

% real_bs = [b0(1) b1(1) b2(1)];
% disp(real_bs)

end