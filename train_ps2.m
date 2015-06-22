clc
clear
global DV IV1 IV2 I J choice_dv se hillbilly

heat = csvread('data5.csv');
DV = heat(:,2);
IV1 = heat(:,3:7)./100;
IV2 = heat(:,8:12)./100;
% Different scaling for different groups: variance in unobserved factors is
% greater for one group of households than for another
hillbilly = heat(:,18); 

I = 900;
J = 5;

choicemat = repmat(DV,1,J);
testmat = repmat(1:J,I,1);
choice_dv = choicemat==testmat;

% b0 = [0 0]
b0 = [0 0 0 0 0 0 0]

% ll(b0)

% options = optimoptions(@fminunc,'Algorithm','quasi-newton')
% optimoptions not supported by Matlab 2012b

% CANNOT USE GRADIENT HERE !!! Gradient of hillbilly is not correct. How to
% calculate ??
options = optimset('LargeScale','on','GradObj','off','Hessian','off')         % LargeScale off is quasi-Newton method in optimset
% options = optimset('Display','iter','LargeScale','off', 'HessUpdate','bfgs' );       % BFGS quasi-Newton

[b, fval,exitflag,output] = fminunc(@train_ps2_fun,b0,options);

disp(['coefficients ' num2str(b) ''])         % cannot print se here ?
disp(['s.e. ' num2str(se') ''])
disp('w: gradient or se of hillbilly is not correct. How to calculate ??')