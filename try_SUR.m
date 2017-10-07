% method see Adv_Econometrics_SysEq_II Slide 11
mydata = csvread('PS0_Data.2.csv',0,0);
mean(mydata)

n = size(mydata,1);
k = 6;
g = 2;

X = zeros(g,k,n);
for i = 1:n
    X(1,1:3,i) = mydata(i,1:3);
    X(2,4:6,i) = mydata(i,4:6);
end

Y = mydata(:,7:8);

% estimate with sols
temp_mat = zeros(k,k);
temp_vec = zeros(k,1);
for j = 1:n
    temp_mat = temp_mat + X(:,:,j)'*X(:,:,j);
    temp_vec = temp_vec + X(:,:,j)'*Y(j,:)';
end
b_sols = inv(1/n.*temp_mat)*(1/n.*temp_vec)

% get residuals of sols
u_sols = zeros(n,g);
for j = 1:n
    u_sols(j,:) = Y(j,:) - (X(:,:,j)*b_sols)';
end
omega_hat = 1/n.*u_sols'*u_sols
inv_omega_hat = inv(omega_hat)

% estimate with fgls
temp_mat2 = zeros(k,k);
temp_vec2 = zeros(k,1);
for j = 1:n
    temp_mat2 = temp_mat2 + X(:,:,j)'*inv_omega_hat*X(:,:,j);
    temp_vec2 = temp_vec2 + X(:,:,j)'*inv_omega_hat*Y(j,:)';
end
b_fgls = inv(1/n.*temp_mat2)*(1/n.*temp_vec2)

% get residuals of fgls
u_fgls = zeros(n,g);
for j = 1:n
    u_fgls(j,:) = Y(j,:) - (X(:,:,j)*b_fgls)';
end

% compute variance-covariance matrix and standard errors
temp_mat3 = zeros(k,k);
for j = 1:n
    temp_mat3 = temp_mat3 + X(:,:,j)'*inv_omega_hat*u_fgls(j,:)'*u_fgls(j,:)*inv_omega_hat*X(:,:,j);
end

Var_fgls = inv(temp_mat2)*temp_mat3*inv(temp_mat2)

std_fgls = diag(Var_fgls)

% w: check kron in answer key
