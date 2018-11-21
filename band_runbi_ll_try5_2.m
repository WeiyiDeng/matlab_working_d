function [b, hessian, grad, standard_error, covariance_matrix, t_stat, exit_flag, output] = band_runbi_ll_try5_2(X, y, beta_0)
global I J dummies se T


IVs = X;

choice_dv = [y 1-y];

beta0 = beta_0;

% options = optimset('LargeScale','on','GradObj','on','Hessian','on','TolFun',1e-6, 'MaxIter',1e4, 'MaxFunEvals', 1e5)         % LargeScale off is quasi-Newton method in optimset
% options = optimset('LargeScale','off','GradObj','off','Hessian','off')
options = optimset('Display','iter','LargeScale','off','GradObj','off','Hessian','off','TolFun',1e-6, 'TolX',1e-14, 'MaxIter',1e3, 'MaxFunEvals', 1e5, 'PlotFcns',@optimplotfirstorderopt);  % ,'OutputFcn', @showJ_history);   %, 'FinDiffType', 'central');

[b, fval,exitflag,output,grad,hessian] = fminunc(@band_bi_ll_try5_2,beta0,options,IVs,choice_dv);

standard_error = sqrt(diag(inv(hessian)));

covariance_matrix = inv(hessian);

t_stat = b./standard_error';

exit_flag = exitflag;

end