function [LL, gr, H] = log_likelihood_gradient_hessian(b)
global IV1 IV2 I J choice_dv se

utility_all = exp(IV1.*b(1)+IV2.*b(2));                              % I*J
pmat = utility_all.*choice_dv./repmat(sum(utility_all,2),1,J);       % I*J    
[r c p] = find(pmat);                                             % I*1
LL = -sum(log(p));                                                % 1*1

format long g
LL

p_all = utility_all./repmat(sum(utility_all,2),1,J);
d = p_all-choice_dv;                     % beware the sequence of what minus what !!

Gt = zeros(I,2);                         % number of variables
Gt(:,1) = sum(IV1.*d, 2);
Gt(:,2) = sum(IV2.*d, 2);
gr = sum(Gt,1);

H = inv(Gt'*Gt)                          % H is actually the inverse of B here (B is approaching -Ht asymptotically) 
disp(gr);
se = sqrt(diag(H));                      % H is the covariance matrix here ? 
disp(se)
fprintf('standard error %d ',se)

end