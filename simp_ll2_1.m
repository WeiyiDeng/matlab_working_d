function LL = simp_ll2_1(beta_0)
global IV I choice_dv

ndraws = 100
b_draws = beta_0(1) + randn(ndraws,1).*beta_0(2);
p_draw = zeros(I,ndraws);
for r = 1:ndraws
    IV = zeros(I,2);
    IV(:,1) = ones(I,1).*b_draws(r);
    utility_all = exp(IV);                   % I*2
    pmat = utility_all.*choice_dv./repmat(sum(utility_all,2),1,2);
    p = pmat(choice_dv ==1);                    % I*1
    p_draw(:,r) = p;                       % 1*1
end
p_avg = mean(p_draw,2);
LL = -sum(log(p_avg));

format long g
LL

end