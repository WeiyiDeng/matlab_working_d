clc
clear all

format long

I = 100000
P = 5

rng(0);
features = randn(I,P-1);
features=[ones(I,1) features-mean(mean(features))];
features(:,3) = features(:,3)-2;

% weeks = randi([-30 30],I,1);
weeks = randi([0 60],I,1);
spaces_empty = randi([1 I],4000,1);
weeks(spaces_empty) = 0;

features(:,5) = weeks;
% norm_weeks = normpdf(weeks,0,8);
norm_weeks = gampdf(weeks,1,8);

% b1 = rand(I,1);
% b2 = rand(I,1).*2+1;

% b1 = randn(1,J).*3
b = [1 2 3 -4 10]'

utilty=[features(:,1:4) norm_weeks]*b;           % utility here is actually the exp of utility 
prob = 1./(1+exp(-utilty));
draw_for_choice=rand(I,1);
choice=draw_for_choice<prob;
% [v,i] = max(prob,[],2);
% choice = i;

DV = [choice 1-choice];
% IV1 = ones(I,J);
IV = features;

b0 = [0 3 0 0 5 3]

options = optimset('disp','iter','LargeScale','on','GradObj','on','Hessian','off','DerivativeCheck','off', 'Jacobian', 'on','TolFun',1e-9, 'MaxIter',1e4, 'MaxFunEvals', 1e6, 'TolX',1e-14)         % LargeScale off is quasi-Newton method in optimset
% options = optimset('LargeScale','on','GradObj','on','Hessian','on')

% options = optimoptions(@fmincon,'Algorithm','interior-point',...
%     'CheckGradients',true,'SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true);

[b, fval,exitflag,output, grad,hessian] = fminunc(@log_likelihood_try_f3,b0,options,DV,IV,I,P);

disp(grad)
disp(full(hessian))
disp(['constants ' num2str(b(1:5)) ''])
disp(['coefficients ' num2str(b(6)) ''])

ww = @(x) 2*exp(-2*x);
plot(gampdf(0.1:0.1:100,1,2))
hold on
plot(ww(0.1:0.1:100))
hold off
