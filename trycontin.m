% obs_mat = [0 0 0; 1 1 0; 0 0 1]
obs_mat = [1 0 0; 0 1 0; 0 0 1];

y = reshape(obs_mat,9,1);

IM = [1 2 3 1 2 3 1 2 3]';
IF = [1 1 1 2 2 2 3 3 3]';

X_mat = [IM IM.^2 IF IF.^2 IM.*IF (IM.*IF).^2];
% IM_quad = IM.^2;
% IF_quad = IF.^2;
% 
% interact = IM.*IF; 

beta0 = [0 0 0 0 0 0 0]';

options = optimset('LargeScale','off','GradObj','off','Hessian','off','TolFun',1e-6, 'TolX',1e-6);

obj_fun = @(b) sum((b(1) + X_mat*b(2:end)-y).^2);

[b,fval] = fminunc(obj_fun, beta0, options)

fitted_y = [ones(length(IM),1) IM IM.^2 IF IF.^2 IM.*IF (IM.*IF).^2]*b
fitted_mat = reshape(fitted_y,3,3)

obs_mat

x1 = 0.01:0.1:3;
x2 = 0.01:0.1:3;
predicted_y = zeros(length(x1),length(x2));
for i = 1:length(x1)
    for j = 1:length(x2)
        predicted_y(i,j) = [1 x1(i) x1(i)^2 x2(j) x2(j)^2 x1(i)*x2(j) (x1(i)*x2(j))^2]*b;
    end
end
surf(predicted_y)        