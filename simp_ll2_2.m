function LL = simp_ll2_2(beta_0)
global IV I choice_dv T

rng(1);
ndraws = 100;
% b_draws = beta_0(1) + randn(ndraws,I).*exp(beta_0(2));
b_draws = beta_0(1) + randn(ndraws,I).*beta_0(2);
draws_panel = repmat(b_draws,1,T);

LL_draw = zeros(ndraws,I*T);
for r = 1:ndraws
    IV = zeros(I*T,2);
    IV(:,1) = ones(I*T,1).*draws_panel(r,:)';
    utility_all = exp(IV);                   % I*2
    pmat = utility_all.*choice_dv./repmat(sum(utility_all,2),1,2);
    p = pmat(choice_dv ==1);                    % I*1
    LL_draw(r,:) = log(p);                       % 1*1
end
time_agg = zeros(ndraws,I);
for i = 1:I
    for t = 1:T
        time_agg(:,i) = time_agg(:,i) + LL_draw(:,(t-1)*I+i);
    end
end
Lid = exp(time_agg);
draw_avg = mean(Lid,1);
LL = -sum(log(draw_avg));

end