function b = use(X)

    function betaHat = ols(data, predIndex, regrIjndices)

x2 = X(:,1);
y = X(:,2);
x1 = ones([size(X,1),1]);
x = [x1,x2];
b = (inv(x'*x))*x'*y;

end

    