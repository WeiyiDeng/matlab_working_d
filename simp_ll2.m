function LL = simp_ll2(beta_0)
global IV I choice_dv

rng(1);
ndraws = 300
b_draws = beta_0(1) + randn(ndraws,I).*beta_0(2);
p_draw = zeros(I,1);
for i = 1:I
    IV = zeros(ndraws,2);
    IV(:,1) = ones(ndraws,1).*b_draws(:,i);
    utility_all = exp(IV);                  
    pmat = utility_all.*repmat(choice_dv(i,:),ndraws,1)./repmat(sum(utility_all,2),1,2);
    [r c p] = find(pmat);                   
    p_draw(i) = mean(p);                    
end
LL = -sum(log(p_draw));

format long g
LL

end