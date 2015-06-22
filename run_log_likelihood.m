clc
clear
global DV IV1 IV2 I J

I = 900;
J = 5;

b0 = [0 0]
heat = csvread('data5.csv');
DV = heat(:,2);
IV1 = heat(:,3:7)./100;
IV2 = heat(:,8:12)./100;

% ll(b0)
choicemat = repmat(DV,1,J);
testmat = repmat(1:J,I,1);
choice_dv = choicemat==testmat;

% options = optimoptions(@fminunc,'Algorithm','quasi-newton')
% optimoptions not supported by Matlab 2012b
options = optimset('LargeScale','off')             % LargeScale off is quasi-Newton method in optimset
% options = optimset('Display','iter','LargeScale','off', 'HessUpdate','bfgs' );       % BFGS quasi-Newton

[b,fval,exitflag,output] = fminunc(@log_likelihood,b0,options);