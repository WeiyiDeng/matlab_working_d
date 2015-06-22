function LL = simp_ll2_3(beta_0)
global IV I choice_dv T

rng(1);
ndraws = 100;
% b_draws = beta_0(1) + randn(ndraws,I).*exp(beta_0(2));
% b_draws = beta_0(1) + randn(ndraws,I).*beta_0(2);
% draws_panel = repmat(b_draws,1,T);

L_draw = zeros(ndraws,I);
for r = 1:ndraws
    idraws = beta_0(1) + randn(1,I).*exp(beta_0(2));
    b_draws = repmat(idraws,T,1);        % T*I
    IV = zeros(T,I,2);
    IV(:,:,1) = ones(T,I).*b_draws;
    utility_all = exp(IV);                   % T*I*2
    pmat = utility_all.*choice_dv./repmat(sum(utility_all,3),[1,1,2]);      %!!!!
    p = sum(pmat,3);                    
    L_draw(r,:) = exp(sum(log(p),1));                       % 1*I
end
avg_draw = mean(L_draw,1);
LL = -sum(log(avg_draw));

end