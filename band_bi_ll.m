function [LL, gr, H] = band_bi_ll(b)
global choice_dv IVs se

% b = (diag([1e-1,1e-1,1e2])*b')';    % w: this does not yield the right result

% const = repmat(b(1),I*T,1);

I = size(IVs,1);
K = size(IVs,2);

const = b(1);

bs = b(2:end)';
FV = IVs*bs;

% exp_util = exp(const+FV);          % utility of choosing the product
% prob=exp_util./(1+exp_util);
exp_util = exp(-(const+FV));         % this is now the utility of the external good
prob=1./(1+exp_util);                % this is still the probability of choosing the product
pmat = [prob 1-prob]; 
pmat = pmat.*choice_dv;
[r c p] = find(pmat);                                             % I*1
LL = -sum(log(p));                                                % 1*1

% format long g
% LL

d = choice_dv(:,1)-prob;                                    % I*1                     
gr = d'*[ones(I,1) IVs];                          % 1*k
gr = -gr;                                                   % w: why reverse ??

H = (repmat(prob.*(1-prob),1,K+1).*[ones(I,1) IVs])'*[ones(I,1) IVs];

% d = prob - choice_dv(:,1);                                % I*1        w: reverse?    
% Gt = repmat(d,1,size(b,2)).*[ones(size(IVs,1),1) IVs];    % I*k
% gr = sum(Gt,1);

% small_value = (2.*rand(1,3)-1)*0.0000005;
% gr = sum(Gt,1).*(1+small_value);

% d = pmat-choice_dv;                                % I*2
% Gt = zeros(size(IVs,1), size(IVs,2)+1);
% Gt(:,1) = sum([zeros(size(IVs,1),1) -ones(size(IVs,1),1)].*d, 2);
% for i = 1:size(IVs,2)
%     Gt(:,i+1) = sum([zeros(size(IVs,1),1) -IVs(:,i)].*d, 2);
% end
% gr = sum(Gt,1);

% % w: hessian not working, why ??
% H = inv(Gt'*Gt);                          % H is actually the inverse of B here (B is approaching -Ht asymptotically)
% disp(gr);
% se = sqrt(diag(H));                      % H is the covariance matrix here ?

%  p_all = utility_all./repmat(sum(utility_all,2),1,J);
%  d = p_all-choice_dv;                     % beware the sequence of what minus what !!
%  
%  Gt = zeros(I,(2*(J-1)+1));                         % number of variables
% %  Gt(:,1) = sum(ones(I,1).*d(:,1), 2);
% %  Gt(:,2) = sum(ones(I,1).*d(:,2), 2);
% %  Gt(:,3) = sum(ones(I,1).*d(:,3), 2);
% %  Gt(:,4) = sum(ones(I,1).*d(:,4), 2);
% for j = 1:(J-1)
%     Gt(:,j) = sum(ones(I,1).*d(:,j), 2);
%     Gt(:,j+1+J-1) = sum(demograph.*d(:,j), 2);
% end
% Gt(:,5) = sum(IV2.*d, 2);
% %  Gt(:,6) = sum(demograph.*d(:,1), 2);
% %  Gt(:,7) = sum(demograph.*d(:,2), 2);
% %  Gt(:,8) = sum(demograph.*d(:,3), 2);
% %  Gt(:,9) = sum(demograph.*d(:,4), 2);
%  
%  gr = sum(Gt,1);
%  
%  H = inv(Gt'*Gt);                          % H is actually the inverse of B here (B is approaching -Ht asymptotically) 
%  disp(gr);
%  se = sqrt(diag(H));                      % H is the covariance matrix here ? 
%  % disp(se);
%  % fprintf('standard error %d ',se);

end