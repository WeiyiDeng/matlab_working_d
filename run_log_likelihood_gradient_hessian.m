clc
clear
global DV IV1 IV2 I J choice_dv se

heat = csvread('data5.csv');
DV = heat(:,2);
IV1 = heat(:,3:7)./100;
IV2 = heat(:,8:12)./100;

I = 900;
J = 5;

choicemat = repmat(DV,1,J);
testmat = repmat(1:J,I,1);
choice_dv = choicemat==testmat;

b0 = [0 0]

% ll(b0)

% options = optimoptions(@fminunc,'Algorithm','quasi-newton')
% optimoptions not supported by Matlab 2012b
options = optimset('LargeScale','on','GradObj','on','Hessian','on')         % LargeScale off is quasi-Newton method in optimset
% options = optimset('Display','iter','LargeScale','off', 'HessUpdate','bfgs' );       % BFGS quasi-Newton

[b, fval,exitflag,output,grad,hessian] = fminunc(@log_likelihood_gradient_hessian,b0,options);

disp(['coefficients ' num2str(b) ''])         % cannot print se here ?
disp(['s.e. ' num2str(se') ''])

% standard_error = sqrt(diag(inv(hessian)))      % hessian provided by fminunc if not supplied