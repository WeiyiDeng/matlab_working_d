% w: to run the estimation
% [x,fval,exitflag,output,grad,hessian] = logregr(imat(:,5), imat(:,7), imat(:,4), 0, 0, 0)

function [x,fval,exitflag,output,grad,hessian] = logregr(l1, l2, hedonic, beta0, beta1, beta2) %Computes likelihood)
bg0=[beta0 beta1 beta2]; % initial values of parameters
options.LargeScale = 'off';
[x,fval,exitflag,output,grad,hessian] = fminunc(@loglike,bg0,options,l1, l2, hedonic)

% calls function loglike