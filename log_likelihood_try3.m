function [LL, gr, H] = log_likelihood_try3(b)
global IV2 I J choice_dv demograph

IV1 = zeros(I,J);
IV1(:,1) = ones(I,1).*b(1);
IV1(:,2) = ones(I,1).*b(2);
IV1(:,3) = ones(I,1).*b(3);
IV1(:,4) = ones(I,1).*b(4);

IV3 = demograph*[b(6:9) 0];

utility_all = exp(IV1+IV2.*b(5)+IV3);                              % I*J
pmat = utility_all.*choice_dv./repmat(sum(utility_all,2),1,J);       % I*J    
[r c p] = find(pmat);                                             % I*1
LL = -sum(log(p));                                                % 1*1

format long g
LL

 p_all = utility_all./repmat(sum(utility_all,2),1,J);
 d = p_all-choice_dv;                     % beware the sequence of what minus what !!
 
 Gt = zeros(I,9);                         % number of variables
 Gt(:,1) = sum(ones(I,1).*d(:,1), 2);
 Gt(:,2) = sum(ones(I,1).*d(:,2), 2);
 Gt(:,3) = sum(ones(I,1).*d(:,3), 2);
 Gt(:,4) = sum(ones(I,1).*d(:,4), 2);
 Gt(:,5) = sum(IV2.*d, 2);
 Gt(:,6) = sum(demograph.*d(:,1), 2);
 Gt(:,7) = sum(demograph.*d(:,2), 2);
 Gt(:,8) = sum(demograph.*d(:,3), 2);
 Gt(:,9) = sum(demograph.*d(:,4), 2);
 
 gr = sum(Gt,1);
 
 H = inv(Gt'*Gt)                          % H is actually the inverse of B here (B is approaching -Ht asymptotically) 
 disp(gr);
 se = sqrt(diag(H));                      % H is the covariance matrix here ? 
 disp(se)
 fprintf('standard error %d ',se)

end