% X is the data
% IV is column index of the predictor variable
% DV is a vector of column indices for the regressor variables
function beta = ols(X,IV,DV)

x2 = X(:,IV);
y = X(:,DV);
x1 = ones([size(X,1),1]);
x = [x1,x2];
beta = (inv(x'*x))*x'*y;

% mylongvariable

end