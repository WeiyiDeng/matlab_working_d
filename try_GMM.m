rng(100)

k = 2;               % # of parameters
n = 50;               
m = 2;               % # of moments

X = 10.*rand(n,k-1);
X = sort(X);
y = 6+X.*3+2.*randn(n,1);
scatter(X,y)

X = [ones(n,1) X];
b_ols = inv(X'*X)*X'*y;
sigma_squared = 1/(n-k)*sum((y-X*b_ols).^2);
Var_b_ols = sigma_squared*inv(X'*X);

b_ols
Var_b_ols

% scatter(X(:,2),y);
% hold on
% x_range = 0:0.1:10;
% plot(x_range,b(1)+x_range.*b(2),'r')
% plot(x_range,(b(1)+2*Var_b(1,1))+x_range.*(b(2)+2*Var_b(2,2)),'--r')
% plot(x_range,(b(1)-2*Var_b(1,1))+x_range.*(b(2)-2*Var_b(2,2)),'--r')
% hold off

% try GMM                                 
% see Applied microeconometrics week1 slide 28
b0_gmm = [1,1]'; 
% g = 1/n.*X'*(y-X*b_gmm);
% W = inv(1/n.*(X'*(y-X*b_gmm)*(X'*(y-X*b_gmm))'));
W0 = eye(2);

options = optimset('Display','iter','LargeScale','off','GradObj','off','Hessian','off','TolFun',1e-6, 'TolX',1e-14, 'MaxIter',1e3, 'MaxFunEvals', 1e5, 'PlotFcns',@optimplotfirstorderopt);  % ,'OutputFcn', @showJ_history);
[b,fval,exitflag,output,grad,hessian] = fminunc(@obj_gmm,b0_gmm,options,X,y,W0,n,k);

b_gmm = b

% computing variance covariance matrix
% w: seems not working properly
g_est = X'*(y-X*b_gmm);
J_est = 1/n.*(g_est*g_est');        % see https://stats.stackexchange.com/questions/112706/two-stage-gmm-estimator-in-matlab
H_est = zeros(m,k);
for j = 1:k
    b_delta = b_gmm;
    b_delta(j) = b_gmm(j)+1.0e-008;
    H_est(:,j) = (1/n.*X'*(y-X*b_delta)-1/n.*g_est)./1.0e-008;
end
V = inv(H_est'*inv(J_est)*H_est)


% ref
% https://stats.stackexchange.com/questions/215604/could-the-covariance-matrix-of-the-moment-conditions-in-gmm-be-ill-conditioned

% https://www.mathworks.com/matlabcentral/fileexchange/12114-gmm?
% - https://www.mathworks.com/matlabcentral/fileexchange/32601-toolkit-on-econometrics-and-economics-teaching?focused=5197665&tab=function

