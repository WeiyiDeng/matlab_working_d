clc
clear all
global DV IV2 I J choice_dv demograph se

I = 10000
J = 5

rng(0);
features = randn(I,J);                                 % eg. costs of transportation (different across individuals and choices)
features=features-mean(mean(features));
% b1 = rand(I,1);
% b2 = rand(I,1).*2+1;
demograph = randn(I,1);                                % eg. income

% b1 = randn(1,J).*3
b1 = [1 2 3 4 5]
b1 = repmat(b1,I,1);
b2 = ones(I,1).*28;
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
% b0 = zeros(1,(2*(J-1)+1))

% options = optimset('LargeScale','on','GradObj','on','Hessian','on','TolFun',1e-6, 'MaxIter',1e4, 'MaxFunEvals', 1e5)         % LargeScale off is quasi-Newton method in optimset
options = optimset('LargeScale','off','GradObj','off','Hessian','off')

[b, fval,exitflag,output,grad,hessian] = fminunc(@log_likelihood_try5,b0,options);

disp(['constants ' num2str(b(1:(J-1))) '']);
disp(['coefficients ' num2str(b(5:(5+J-1))) '']);
% disp(['standard errors ' num2str(se') '']);
% t_stat = b(5:(5+J-1))./se(5:(5+J-1))';
% disp(['t statistics ' num2str(t_stat) '']);

standard_error = sqrt(diag(inv(hessian)))      
% These s.e. are correct only when grad and hessian are NOT provided by me 
% (turn off GradObj and Hessian to use these standard errors)
% w: why is the case ??

t_stat = b(5:(5+J-1))./standard_error(5:(5+J-1))';
disp(['t statistics ' num2str(t_stat) '']);
