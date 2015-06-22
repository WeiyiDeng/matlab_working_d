clc
clear
global I choice_dv T

I = 1000
% b0 = [5 exp(2)]
b0 = [0 0]
% b0 = 5
T = 10

rng(0);
% beta0 = 0.6
% mu_j = repmat(beta0,I,1)

beta0 = 12; 
beta1 = 1;

mu_j = beta0 + randn(I,1).*beta1;
mu_jmat = repmat(mu_j,T,1);
mu_k = zeros(I,1)
mu_kmat = repmat(mu_k,T,1);
muj = reshape(mu_jmat,T*I,1);
muk = reshape(mu_kmat,T*I,1);
prob = exp(muj)./(exp(muj)+exp(muk));         % I*1
draw_for_choice=rand(I*T,1);
choice=prob>draw_for_choice;

choice_dv = [choice, 1-choice];                % I*2

options = optimset('LargeScale','off','GradObj','off','Hessian','off','display','iter')

[beta_0, fval,exitflag,output,grad,hessian] = fminunc(@simp_ll2_2,b0,options);

se = sqrt(diag(inv(hessian))) 
