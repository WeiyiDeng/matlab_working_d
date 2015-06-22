clc
clear
global I choice_dv

I = 10000
b0 = 0

rng(0);
beta0 = 5
mu_j = repmat(beta0,I,1)
mu_k = zeros(I,1)
prob = exp(mu_j)./(exp(mu_j)+exp(mu_k))         % I*1
draw_for_choice=rand(I,1)
choice=prob>draw_for_choice

choice_dv = [choice, 1-choice]                % I*2

options = optimset('LargeScale','off','GradObj','off','Hessian','off')

[beta_0, fval,exitflag,output,grad,hessian] = fminunc(@simp_ll,b0,options);

se = sqrt(diag(inv(hessian))) 
