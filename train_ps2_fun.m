function [LL, gr, H] = train_ps2_fun(b)
global IV1 IV2 I J choice_dv se hillbilly

m = (b(3).*hillbilly)+(1-hillbilly)
m = repmat(m,1,J)
CONST = zeros(I,J);
CONST(:,1) = ones(I,1).*b(4);
CONST(:,2) = ones(I,1).*b(5);
CONST(:,3) = ones(I,1).*b(6);
CONST(:,4) = ones(I,1).*b(7);

utility_all = exp(CONST+IV1.*b(1)+IV2.*b(2)).^m
% utility_all = exp(IV1.*b(1)+IV2.*b(2));                              % I*J
pmat = utility_all.*choice_dv./repmat(sum(utility_all,2),1,J);       % I*J    
[r c p] = find(pmat);                                             % I*1
LL = -sum(log(p));                                                % 1*1

format long g
LL

p_all = utility_all./repmat(sum(utility_all,2),1,J);
d = p_all-choice_dv;                     % beware the sequence of what minus what !!

Gt = zeros(I,7);                         % number of variables
Gt(:,1) = sum(IV1.*d, 2);
Gt(:,2) = sum(IV2.*d, 2);
Gt(:,3) = sum(hillbilly.*d(choice_dv==1), 2);          % gradient for hillbilly is not correct. How to calculate ??
Gt(:,4) = sum(ones(I,1).*d(:,1), 2);
Gt(:,5) = sum(ones(I,1).*d(:,2), 2);
Gt(:,6) = sum(ones(I,1).*d(:,3), 2);
Gt(:,7) = sum(ones(I,1).*d(:,4), 2);
gr = sum(Gt,1);

H = inv(Gt'*Gt)                          % H is actually the inverse of B here (B is approaching -Ht asymptotically) 
disp(gr);
se = sqrt(diag(H));                      % H is the covariance matrix here ? 
disp(se)
fprintf('standard error %d ',se)

end