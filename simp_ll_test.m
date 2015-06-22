function LL = simp_ll_test(beta_0)
global IV I choice_dv T

ndraws = 100;
% b_draws = beta_0(1) + randn(ndraws,I).*exp(beta_0(2));
% b_draws = beta_0(1) + randn(ndraws,I).*beta_0(2);
% draws_panel = repmat(b_draws,1,T);

threshold = 4;
% c = 1;
L_draw = zeros(ndraws,I);
seed = 1;
% rng(seed);
for r = 1:ndraws
    % Accept-reject draws from normal
    rng(seed);
    trun_randn = randn(1,I);
    if sum(abs(trun_randn)>threshold)>=1
        trun_randn(abs(trun_randn)>threshold) = [];
    end
    while length(trun_randn)<I
        seed = seed + 1;                 % w: necessary ????????????/
        trun_randn = [trun_randn randn(1,I-length(trun_randn))];
        if sum(abs(trun_randn)>threshold)>=1
            trun_randn(abs(trun_randn)>threshold) = [];
        end
    end
    % idraws = beta_0(1) + randn(1,I).*exp(beta_0(2));
    idraws = beta_0(1) + trun_randn.*exp(beta_0(2));         % 1*I
    % idraws = beta_0(1) + randn(1,I).*(c./(1+exp(beta_0(2))));
    b_draws = repmat(idraws,T,1);        % T*I
    IV = zeros(T,I,2);
    IV(:,:,1) = ones(T,I).*b_draws;
    utility_all = exp(IV);                   % T*I*2
    pmat = utility_all.*choice_dv./repmat(sum(utility_all,3),[1,1,2]);    
    p = sum(pmat,3);                    
    L_draw(r,:) = exp(sum(log(p),1));                       % 1*I
    seed = seed+1;
end
avg_draw = mean(L_draw,1);
LL = -sum(log(avg_draw));

end