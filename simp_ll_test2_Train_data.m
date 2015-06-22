function LL = simp_ll_test2_Train_data(beta_0)
global I IV choice_dv J N K Respondent_mat

D = 100;

% bb = beta_0(2*J-1:2*(J-1+K));
% b = reshape(bb,2,K);                            % K*2
b = reshape(beta_0,2,K)';

mu = b(:,1);
se = b(:,2);

seed = 1;
prob_simu = zeros(I,1);
Respondent = Respondent_mat(:,1);

% ncount = 0;                                  % number of observations above
for i = 1:I
    rng(seed);
    idraws = repmat(mu,1,D) + randn(K,D).*repmat(se,1,D);        % K*D        nVar*nDraws
    [r c p] = find(Respondent==i);
    iIV = IV(r,:,:);                               % T*J*K
    iTime = length(r);
    ichoices = choice_dv(r,:);             % T*J
    % p_seq = ones(1,D);
    sum_pt = zeros(1,D);
    for t = 1:iTime
        expU_t = exp(squeeze(iIV(t,:,:))*idraws);            % J*D
        p_chosen_t = expU_t(ichoices(t,:)==1,:)./sum(expU_t,1);
        sum_pt = sum_pt + log(p_chosen_t);
        % p_seq = p_seq.*p_chosen_t;
    end
    p_seq = exp(sum_pt);
    prob_simu(i,1) = mean(p_seq);
    seed = seed+1;
end

LL = -sum(log(prob_simu));

end