% simulation Lec 11 Slide 45-46 Batch means method example
% Mar 15
% Apr 4 replication

X = [2.69
2.08
3.63
7.83
2.36
1.64
4.18
5.84
4.85];

mean_X = mean(X)
SE_X = sqrt(sum((A - mean(A)).^2)/(length(A)-1))
SE_X_bar = SE_X/sqrt(length(A))          % var(X_bar) = var(Xi)/n

CI_X_bar_low = mean_X-2*SE_X_bar
CI_X_bar_high = mean_X+2*SE_X_bar