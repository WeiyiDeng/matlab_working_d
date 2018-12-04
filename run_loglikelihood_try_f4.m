clc
clear all

format long

I = 500000
P = 5

rng(0);
features = randn(I,P-1);
features=[ones(I,1) features-mean(mean(features))];
features(:,3) = features(:,3)-2;

weeks = randi([-30 30],I,1);
spaces_empty = randi([1 I],4000,1);
weeks(spaces_empty) = 0;

features(:,5) = weeks;
norm_weeks = normpdf(weeks,0,8);

diff_weeks = randi([1 423],I,1);
spaces_empty_diff = randi([1 I],3000,1);
% diff_weeks(spaces_empty_diff) = 0;

features(:,4) = diff_weeks;
norm_weeks_diff = gampdf(diff_weeks,1.25,6);           %w% gamma_derivative_wrt_alpha not defined for x=0(returns NAN)
% norm_weeks_diff = normpdf(diff_weeks,0,60);       
% sth = normpdf(diff_weeks,0,6);

% b1 = rand(I,1);
% b2 = rand(I,1).*2+1;

% b1 = randn(1,J).*3
b = [1 -2 3 -15 10]'

utilty=[features(:,1:3) norm_weeks_diff norm_weeks]*b;           % utility here is actually the exp of utility 
prob = 1./(1+exp(-utilty));
draw_for_choice=rand(I,1);
choice=draw_for_choice<prob;
% [v,i] = max(prob,[],2);
% choice = i;

DV = [choice 1-choice];
% IV1 = ones(I,J);
IV = features;

b0 = [0 -3 5 -10 8 2 0 2]

options = optimset('disp','iter','LargeScale','on','GradObj','on','Hessian','off','DerivativeCheck','off', 'Jacobian', 'on','TolFun',1e-9, 'MaxIter',3e2, 'MaxFunEvals', 1e6, 'TolX',1e-14)         % LargeScale off is quasi-Newton method in optimset
% options = optimset('LargeScale','on','GradObj','on','Hessian','on')

% options = optimoptions(@fmincon,'Algorithm','interior-point',...
%     'CheckGradients',true,'SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true);

[b, fval,exitflag,output, grad,hessian] = fminunc(@log_likelihood_try_f4,b0,options,DV,IV,I,P);

disp(grad)
disp(full(hessian))
% disp(['constants ' num2str(b(1:5)) ''])
% disp(['coefficients ' num2str(b(6)) ''])
disp(b)

plot(gampdf(0.1:0.1:100,1,6))
hold on
plot(normpdf(0.1:0.1:100,0,8))
hold off
