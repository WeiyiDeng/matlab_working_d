clc
clear all
global DV IV2 I J choice_dv demograph

I = 1000
J = 5

rng(0);
features = randn(I,J);
features=features-mean(mean(features));
% b1 = rand(I,1);
% b2 = rand(I,1).*2+1;
demograph = randn(I,1);

% b1 = randn(1,J).*3
b1 = [1 2 3 4 5]
b1 = repmat(b1,I,1);
b2 = ones(I,1).*10;
b2 = repmat(b2,1,J);
b3 = [1 3 2 4 6]

utilty=exp(b1+ b2.*features + demograph*b3);           % utility here is actually the exp of utility 
utilty_other=sum(utilty,2)
utilty_other=repmat(utilty_other,1,J);
prob=utilty./utilty_other;
prob=cumsum(prob')';
draw_for_choice=rand(I,1);
draw_for_choice=repmat(draw_for_choice,1,J);
choice=prob<draw_for_choice;
choice=sum(choice,2)+1;
% [v,i] = max(prob,[],2);
% choice = i;

DV = choice;
choicemat = repmat(DV,1,J);
testmat = repmat(1:J,I,1);
choice_dv = choicemat==testmat;
% IV1 = ones(I,J);
IV2 = features;

b0 = [0 0 0 0 0 0 0 0 0]

% options = optimset('LargeScale','on','GradObj','on','Hessian','on','TolFun',1e-6, 'MaxIter',1e4, 'MaxFunEvals', 1e5)         % LargeScale off is quasi-Newton method in optimset
options = optimset('LargeScale','on','GradObj','on','Hessian','on')

[b, fval,exitflag,output] = fminunc(@log_likelihood_try3,b0,options);

disp(['constants ' num2str(b(1:4)) ''])
disp(['coefficients ' num2str(b(5:9)) ''])




