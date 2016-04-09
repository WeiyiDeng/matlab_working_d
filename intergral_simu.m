% % Lec05 E2
% % take integral of (x^2+x) from 1 to 3
% syms x
% x = 3;
% intp1 = x^3/3+x^2/2;
% x = 1;
% intp2 = x^3/3+x^2/2;
% analytical_int = intp1-intp2

k = 10000;
d = rand(k,1).*2+1;
g = d.^2+d;
theta = mean(g)*2

% Lec05 E3
k = 1000;
x = rand(k,1);
y = rand(k,1);
g = 4*x.^2.*y+y.^2;
theta = mean(g)

% Calculus textbook Po310 E2
 % take integral of (x^2*y^2) in (0<x<2, -1<y<1)
k = 1000;
x = rand(k,1).*2;
y = rand(k,1).*2-1;
g = x.^2.*y.^2;
theta = 4*mean(g)

% Calculus textbook Po313 E7
k = 1000;
x = rand(k,1);
y = rand(k,1);
ind = zeros(k,1);
ind(y.^2>x) = 1;
D = mean(ind)*1
x_star = x(y.^2>x);
y_star = y(y.^2>x);
g = exp(x_star./y_star);
theta = D*mean(g)