clc
clear
global I choice_dv T

I = 200
% b0 = [5 exp(2)]
b0 = [0 exp(0)]
% b0 = 5
T = 50

rng(0);
% beta0 = 0.6
% mu_j = repmat(beta0,I,1)

beta0 = 7; 
beta1 = 1;

mu_j = beta0 + randn(I,1).*beta1;            % I*1
muj = repmat(mu_j,T,1);                  % IT*1
mu_k = zeros(I,1);
muk = repmat(mu_k,T,1);
% muj = reshape(mu_jmat,T*I,1);
% muk = reshape(mu_kmat,T*I,1);
prob = exp(muj)./(exp(muj)+exp(muk));         
draw_for_choice=rand(I*T,1);
choice=prob>draw_for_choice;

choice_dv = zeros(T,I,2);                     % number of alternatives
choice_dv(:,:,1) = reshape(choice,T,I);
choice_dv(:,:,2) = 1-reshape(choice,T,I);

% choice_dv = [choice, 1-choice];                % I*2

options = optimset('LargeScale','off','GradObj','off','Hessian','off','display','iter')

[beta_0, fval,exitflag,output,grad,hessian] = fminunc(@simp_ll2_3,b0,options);

se = sqrt(diag(inv(hessian))) 
