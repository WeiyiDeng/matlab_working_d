X = 10.*rand(50,1);
X = sort(X);
y = 6+X.*3+10.*randn(50,1);
scatter(X,y)

X = [ones(50,1) X];
b = inv(X'*X)*X'*y;

k = size(X,2)-1;
n = size(X,1);
sigma_squared = 1/(n-k)*sum((y-X*b).^2);
Var_b = sigma_squared*inv(X'*X);

b
Var_b

scatter(X(:,2),y);
hold on
x_range = 0:0.1:10;
plot(x_range,b(1)+x_range.*b(2),'r')
plot(x_range,(b(1)+2*Var_b(1,1))+x_range.*(b(2)+2*Var_b(2,2)),'--r')
plot(x_range,(b(1)-2*Var_b(1,1))+x_range.*(b(2)-2*Var_b(2,2)),'--r')
hold off