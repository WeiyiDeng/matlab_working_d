% function [b, hessian, grad, standard_error, covariance_matrix, t_stat, exit_flag, output] = band_runbi_ll_p(X, y, dummy_X, beta_0)
function [b, hessian, grad, standard_error, covariance_matrix, t_stat, exit_flag, output] = band_runbi_ll_trend_format(X, y, beta_0, innov_X, explor_X)
global I J dummies se T

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
choice_dv = [y 1-y];

% plus = zeros(length(choice_dv),1);
% plus(IVs(:,2)==0)=100000;
% IVs(:,2) = IVs(:,2)+plus;
% clearvars plus

% week_IV = gampdf(IVs(:,2),0.9337,27.4738);              % w: NOTICE new lines here for fixed gamma parameters
% week_IV(IVs(:,2)<1)=0;      
% 
% innov_WD_multip = zeros(size(innov_X));
% for i = 1:size(innov_X,2);
%     innov_WD_multip(:,i) = innov_X(:,i).*week_IV;
% end
% 
% explor_WD_multip = zeros(size(explor_X));
% for j = 1:size(explor_X,2);
%     explor_WD_multip(:,j) = explor_X(:,j).*week_IV;
% end

% beta0 = [0 0 0]
beta0 = beta_0;

% options = optimset('LargeScale','on','GradObj','on','Hessian','on','TolFun',1e-6, 'MaxIter',1e4, 'MaxFunEvals', 1e5)         % LargeScale off is quasi-Newton method in optimset
% options = optimset('LargeScale','off','GradObj','off','Hessian','off')
options = optimset('Display','iter','LargeScale','off','GradObj','off','Hessian','off','TolFun',1e-6, 'TolX',1e-14, 'MaxIter',1e3, 'MaxFunEvals', 1e5, 'PlotFcns',@optimplotfirstorderopt);  % ,'OutputFcn', @showJ_history);   %, 'FinDiffType', 'central');

% [b, fval,exitflag,output,grad,hessian] = fminunc(@band_bi_ll_i2,beta0,options,IVs,choice_dv, innov_X, explor_X, week_IV, innov_WD_multip, explor_WD_multip);
[b, fval,exitflag,output,grad,hessian] = fminunc(@band_bi_ll_trend_format,beta0,options,IVs,choice_dv, innov_X, explor_X);

% disp(['constant ' num2str(b(1)) '']);
% disp(['coefficients ' num2str(b(2:end)) '']);

standard_error = sqrt(diag(inv(hessian)));
% These s.e. are correct only when grad and hessian are NOT provided by me 
% (turn off GradObj and Hessian to use these standard errors)
% w: why is the case ??

covariance_matrix = inv(hessian);

t_stat = b./standard_error';
% disp(['t statistics ' num2str(t_stat) '']);

exit_flag = exitflag;

% real_bs = [b0(1) b1(1) b2(1)];
% disp(real_bs)

end