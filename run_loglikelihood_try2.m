clc
clear
global DV I J choice_dv

I = 10000
J = 5

rng(0);
features = randn(I,J);
% b1 = rand(I,1);
% b2 = rand(I,1).*2+1;

% b1 = randn(1,J).*3
b1 = [1 2 3 4 5]
b1 = repmat(b1,I,1);
% b2 = ones(I,1).*3;
% b2 = repmat(b2,1,J);

% utilty=exp(b1+ b2.*features);           % utility here is actually the exp of utility 
utilty=exp(b1);
utilty_other=sum(utilty,2)
utilty_other=repmat(utilty_other,1,J);
prob=utilty./utilty_other;
prob=cumsum(prob')';
draw_for_choice=rand(I,1);
draw_for_choice=repmat(draw_for_choice,1,J);
choice=prob<draw_for_choice;
choice=sum(choice,2)+1;

DV = choice;
choicemat = repmat(DV,1,J);
testmat = repmat(1:J,I,1);
choice_dv = choicemat==testmat;
% IV1 = ones(I,J);
% IV2 = features;

b0 = [0 0 0 0]

options = optimset('LargeScale','on','GradObj','off','Hessian','off')         % LargeScale off is quasi-Newton method in optimset

[b, fval,exitflag,output] = fminunc(@log_likelihood_try2,b0,options);

disp(['constants ' num2str(b(1:4)) ''])
% disp(['coefficients ' num2str(b(5)) ''])




