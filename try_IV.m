mydata = csvread('PS0_Data.1.csv',0,0);
mean(mydata)

X = mydata(:,1:3);
y1 = mydata(:,4);
y2 = mydata(:,5);
z1 = mydata(:,6);
z2 = mydata(:,7);

% exact;y identified case
X_step1 = [X z1];
b_step1 = inv(X_step1'*X_step1)*X_step1'*y2
y2_predict = X_step1*b_step1;

X_step2 = [X y2_predict];
b_step2 = inv(X_step2'*X_step2)*X_step2'*y1

n = size(X,1);
k = 4;

XX = [X y2];
std = sqrt(diag(mean((y1-XX*b_step2).^2)*((X_step2'*XX)\eye(size(XX,2)))))          % w: check?


% over identified case
X_step1 = [X z1 z2];
b_step1 = inv(X_step1'*X_step1)*X_step1'*y2
y2_predict = X_step1*b_step1;

X_step2 = [X y2_predict];
b_step2 = inv(X_step2'*X_step2)*X_step2'*y1

n = size(X,1);
k = 4;

XX = [X y2];
std = sqrt(diag(mean((y1-XX*b_step2).^2)*((X_step2'*XX)\eye(size(XX,2)))))          % w: check?

%%
% try GMM over identified case
m = 5;                      % # of moment conditions

b0_gmm = ones(k,1);
W0 = eye(m);
X_gmm = [X y2];
y_gmm = y1;
Z = [X z1 z2];

options = optimset('Display','iter','LargeScale','off','GradObj','off','Hessian','off','TolFun',1e-6, 'TolX',1e-14, 'MaxIter',1e3, 'MaxFunEvals', 1e5, 'PlotFcns',@optimplotfirstorderopt);  % ,'OutputFcn', @showJ_history);
[b,fval,exitflag,output,grad,hessian] = fminunc(@obj_gmm_IV,b0_gmm,options,X_gmm,y_gmm,W0,Z,n,k);

b0_gmm = b
W = inv(1/n.*(Z'*(y_gmm-X_gmm*b0_gmm)*(Z'*(y_gmm-X_gmm*b0_gmm))'))             % W computation incorrect? Error Message: Matrix is close to singular or badly scaled

W0 = W;
[b,fval,exitflag,output,grad,hessian] = fminunc(@obj_gmm_IV,b0_gmm,options,X_gmm,y_gmm,W0,Z,n,k);
b
