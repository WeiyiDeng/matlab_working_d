clc
clear
global I IV DV Time J N K

xdata = csvread('xdata.csv');
ydata = csvread('ydata.csv');
Time = csvread('times.csv');

X = zeros(4308,28);            % number of observations * 7 variables for the 4 alternatives
r = 0;
while r < 4308
    X(r+1,:) = reshape(xdata((7*r+1:7*r+7),:)',1,[]);
    r = r+1
end

DV = ydata;

N = 4308             % number of observations
J = 4                % number of alternatives
I = 361              % number of individuals
K = 6                % number of variables

IV1 = reshape(X(:,1:4)',N*J,1);            % vectorize: first observation with 4 alternatives, second observation with 4 alternatives...
IV2 = reshape(X(:,5:8)',N*J,1);
IV3 = reshape(X(:,9:12)',N*J,1);
IV4 = reshape(X(:,13:16)',N*J,1);
IV5 = reshape(X(:,21:24)',N*J,1);
IV6 = reshape(X(:,25:28)',N*J,1);
IV = [IV1 IV2 IV3 IV4 IV5 IV6];

tic

% b0 = zeros(1,12)
b0 = [-.857 0.4 -.183 0.1 2.098 1 1.525 0.8 -8.285 2 -8.53 2];

% options = optimset('LargeScale','off','GradObj','off','Hessian','off','display','iter')
options = optimset('LargeScale','off','GradObj','off','Hessian','off','display','iter', 'MaxIter',1e4, 'MaxFunEvals', 1e5)         % LargeScale off is quasi-Newton method in optimset

[beta_0, fval,exitflag,output,grad,hessian] = fminunc(@try_mxl,b0,options);

beta_0
se = sqrt(diag(inv(hessian)))
t_stat = beta_0./se'

toc