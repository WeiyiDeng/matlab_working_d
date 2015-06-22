clc
clear

b0 = [0 0]
% I = 900
% J = 5
% ll(b0)

% options = optimoptions(@fminunc,'Algorithm','quasi-newton')
% optimoptions not supported by Matlab 2012b
options = optimset('LargeScale','off')             % LargeScale off is quasi-Newton method in optimset
% options = optimset('Display','iter','LargeScale','off', 'HessUpdate','bfgs' );       % BFGS quasi-Newton

[b,fval,exitflag,output] = fminunc(@ll,b0,options);
