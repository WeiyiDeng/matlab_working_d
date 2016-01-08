X = 10.*rand(50,1);
X = sort(X);
y = 6+X.*3+10.*randn(50,1);
scatter(X,y)

X = [ones(50,1) X];
b = inv(X'*X)*X'*y;

scatter(X(:,2),y);
hold on
x_range = 0:0.1:10;
plot(x_range,b(1)+x_range.*b(2),'r')
hold off