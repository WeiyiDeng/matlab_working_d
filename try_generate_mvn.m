x = randn(500,1);
epsilon = randn(500,1);
y = x+epsilon;
b = inv(x'*x)*x'*y
e = y - b.*x;
e_hat = mean(e)
s_square = 1/(500-2)*sum((e-e_hat).^2)

% try mvn with cov between variables not equal to 0
mu = [0 0]
sigma = [1 0.5;0.5 1]           % cov = 0.5
% sigma = [1 0;0 1]             % cov = 0
trymvn = mvnrnd(mu, sigma, 500);
x = trymvn(:,1);
epsilon = trymvn(:,2);
y = x+epsilon;
b = inv(x'*x)*x'*y
e = y - b.*x;
e_hat = mean(e)
s_square = 1/(500-2)*sum((e-e_hat).^2)