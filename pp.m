% Weiyi Deng 401716
load food.dat
X = food;
B = use(X)

syms m
Y = B'*[1;m]
scatter(X(:,1),X(:,2),'black')
hold all
h = ezplot(Y,[0,max(X(:,1))]);
set(h,'Color','red')
xlabel('weekly income')
ylabel('food expenditure')
text(max(X(:,1))*0.6,max(X(:,2))*0.6,['y=',num2str(B(1,1)),'+',num2str(B(2,1)),'*x'],'Color','blue')
title('OLS results')
hold off