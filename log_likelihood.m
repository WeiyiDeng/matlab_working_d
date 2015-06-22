function LL = log_likelihood(b)
global DV IV1 IV2 I J choice_dv

utility_all = exp(IV1.*b(1)+IV2.*b(2));
pmat = utility_all.*choice_dv./repmat(sum(utility_all,2),1,J);
[r c p] = find(pmat);
LL = -sum(log(p));

format long g
LL

end