function LL = try_mxl(beta_0)
global I IV DV Time J N K

D = 100;                                     % number of draws
prob_simu = zeros(I,1);
b = reshape(beta_0,2,K)';                    % K*2
mu = b(:,1);
se = b(:,2);

seed = 1;
% rng(seed);
ncount = 0;                                  % number of observations above
for i = 1:I
    rng(seed);
    idraws = repmat(mu,1,D) + randn(K,D).*repmat(se,1,D);        % K*D        nVar*nDraws
    iIV = IV(ncount*J+1:ncount*J+Time(i)*J,:);                  %TJ*K
    exp_utility = exp(iIV*idraws);
    p_seq = ones(1,D);
    for t = 1:Time(i)
        expU_t = exp_utility((t-1)*J+1:t*J,:);
        p_chosen_t = expU_t(DV(ncount+t),:)./sum(expU_t,1);
        p_seq = p_seq.*p_chosen_t;
    end
    prob_simu(i,1) = mean(p_seq);
    ncount = ncount+Time(i);
    seed = seed+1;
end

LL = -sum(log(prob_simu));

end
        
        
        