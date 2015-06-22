function [LL gr] = log_likelihood_gradient(b)
global IV1 IV2 I J choice_dv

utility_all = exp(IV1.*b(1)+IV2.*b(2));
pmat = utility_all.*choice_dv./repmat(sum(utility_all,2),1,J);
[r c p] = find(pmat);
LL = -sum(log(p));

format long g
LL

p_all = utility_all./repmat(sum(utility_all,2),1,J);
d = p_all-choice_dv;                     % beware the sequence of what minus what !!

Gt = zeros(I,2);                        % number of variables
Gt(:,1) = sum(IV1.*d, 2);
Gt(:,2) = sum(IV2.*d, 2);
gr = sum(Gt,1);

end