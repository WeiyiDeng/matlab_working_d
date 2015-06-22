clc
clear
global DV IV1 IV2 I J choice_dv

b0 = [0 0]
heat = csvread('data5.csv');
DV = heat(:,2);
IV1 = heat(:,3:7)./100;
IV2 = heat(:,8:12)./100;

I = 900;
J = 5;

choicemat = repmat(DV,1,J);
testmat = repmat(1:J,I,1);
choice_dv = choicemat==testmat;

% ll(b0)

% options = optimoptions(@fminunc,'Algorithm','quasi-newton')
% optimoptions not supported by Matlab 2012b
options = optimset('LargeScale','off','GradObj','on')         % LargeScale off is quasi-Newton method in optimset
% options = optimset('Display','iter','LargeScale','off', 'HessUpdate','bfgs' );       % BFGS quasi-Newton

[b,fval,exitflag,output] = fminunc(@log_likelihood_gradient,b0,options);