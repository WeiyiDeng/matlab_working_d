function [LL, g, H] = log_likelihood_try_f5(b,DV,IV,I,P)

b_basics = b(1:5);

exp_util = exp(-[IV(:,1:3) gampdf(IV(:,4),exp(b(7)),exp(b(8))) normpdf(IV(:,5),0,exp(b(6)))]*b_basics');         % this is now the utility of the external good
% exp_util = exp(-[IV(:,1:3) normpdf(IV(:,4),0,exp(b(7))) normpdf(IV(:,5),0,exp(b(6)))]*b_basics'); 
prob=1./(1+exp_util);                % this is still the probability of choosing the product
pmat = [prob 1-prob]; 
pmat = pmat.*DV;
[r c p] = find(pmat);                                             % I*1
LL = -sum(log(p));                                          % 1*1

g(1:3) = IV(:,1:3)'*(prob-DV(:,1));
g(4) = gampdf(IV(:,4),exp(b(7)),exp(b(8)))'*(prob-DV(:,1));
% g(4) = normpdf(IV(:,4),1,exp(b(7)))'*(prob-DV(:,1));
g(5) = normpdf(IV(:,5),0,exp(b(6)))'*(prob-DV(:,1));
% normal_derivative_wrt_sigmas = b(end-1).*2^(-5/2).*(IV(:,5).^2-1).*normpdf(IV(:,5),0,exp(b(end)))./exp(b(end));
normal_derivative_wrt_sigmas = normal_derivative(IV(:,5), 0, b(6)); 
g(6) = b(5).*normal_derivative_wrt_sigmas'*(prob-DV(:,1));
gamma_derivative_alpha = gamma_derivative_wrt_alpha(IV(:,4), b(7), b(8)); 
g(7) = b(4).*gamma_derivative_alpha'*(prob-DV(:,1));
gamma_derivative_beta = gamma_derivative_wrt_beta(IV(:,4), b(7), b(8)); 
g(8) = b(4).*gamma_derivative_beta'*(prob-DV(:,1));
% gamma_derivative_beta = normal_derivative(IV(:,4), 0, b(7)); 
% g(7) = b(4).*gamma_derivative_beta'*(prob-DV(:,1));

% H = IV'*diag(prob*(1-prob)')*IV;

%  p_all = utility_all./repmat(sum(utility_all,2),1,J);
%  d = p_all-choice_dv;                     % beware the sequence of what minus what !!
%  
%  Gt = zeros(I,5);                         % number of variables
%  Gt(:,1) = sum(ones(I,1).*d(:,1), 2);
%  Gt(:,2) = sum(ones(I,1).*d(:,2), 2);
%  Gt(:,3) = sum(ones(I,1).*d(:,3), 2);
%  Gt(:,4) = sum(ones(I,1).*d(:,4), 2);
%  Gt(:,5) = sum(IV2.*d, 2);
%  gr = sum(Gt,1);
%  
%  H = inv(Gt'*Gt)                          % H is actually the inverse of B here (B is approaching -Ht asymptotically) 
%  disp(gr);
%  se = sqrt(diag(H));                      % H is the covariance matrix here ? 
%  disp(se)
%  fprintf('standard error %d ',se)

end