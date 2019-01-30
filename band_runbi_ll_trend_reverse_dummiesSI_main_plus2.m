% function [b, hessian, grad, standard_error, covariance_matrix, t_stat, exit_flag, output] = band_runbi_ll_p(X, y, dummy_X, beta_0)
function [b, hessian, grad, standard_error, covariance_matrix, t_stat, exit_flag, output] = band_runbi_ll_trend_reverse_dummiesSI_main_plus2(X, trend_hat, band_age, topics_count, dummy_prep, None0s_X_N, S, y, beta_0)
global I J dummies se T

IVs = X;

choice_dv = [y 1-y];

beta0 = beta_0;

% options = optimset('LargeScale','on','GradObj','on','Hessian','on','TolFun',1e-6, 'MaxIter',1e4, 'MaxFunEvals', 1e5)         % LargeScale off is quasi-Newton method in optimset
% options = optimset('LargeScale','off','GradObj','off','Hessian','off')
options = optimset('Display','iter','LargeScale','off','GradObj','off','Hessian','off','TolFun',1e-6, 'TolX',1e-14, 'MaxIter',1e3, 'MaxFunEvals', 1e5, 'PlotFcns',@optimplotfirstorderopt);  % ,'OutputFcn', @showJ_history);   %, 'FinDiffType', 'central');

% [b, fval,exitflag,output,grad,hessian] = fminunc(@band_bi_ll_i2,beta0,options,IVs,choice_dv, innov_X, explor_X, week_IV, innov_WD_multip, explor_WD_multip);
[b, fval,exitflag,output,grad,hessian] = fminunc(@band_bi_ll_trend_reverse_dummiesSI_main_plus2,beta0,options,IVs, trend_hat, band_age, topics_count, dummy_prep, None0s_X_N, S, choice_dv);

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