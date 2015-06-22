syms x
y = exp(-x^2/2)
f = int(y, x, -inf, inf)   % 求积分integral
f              % 得正态分布pdf里的常数