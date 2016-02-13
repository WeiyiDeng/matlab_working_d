% w: matlab cov() gives different result from calculated as in AdvStat course, why??
X = [9 1;5 3;1 2]
mu = mean(X)
transmu = repmat(mu,3,1)
varcov = 1/3.*[(X- transmu)'*(X- transmu)]
x1 = [9 1 5]
x2 = [1 3 2]
cov(x1,x2)

% if cov() is using 1/(n-1) the following result is still not the same
1/2.*[(X- transmu)'*(X- transmu)]
var(x1)
var(x2)