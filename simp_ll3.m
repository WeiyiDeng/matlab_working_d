function LL = simp_ll3(beta_0)
global IV I choice_dv

IV = zeros(I,2);
IV(:,1) = ones(I,1).*beta_0;
utility_all = exp(IV);                   % I*2
pmat = utility_all.*choice_dv./repmat(sum(utility_all,2),1,2);
p = pmat(choice_dv ==1);                    % I*1
LL = -sum(log(p));                       % 1*1

format long g
LL

end