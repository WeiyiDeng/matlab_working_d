function LL = mnl_ll(beta_0)
global I IV choice_dv J K Respondent_mat

D = 100;                                        % number of draws

bb = beta_0(2*J-1:2*(J-1+K));
b = reshape(bb,2,K)';                            % K*2

mu = b(:,1);
se = b(:,2);

bc = beta_0(1:2*J-2);
c = reshape(bc,2,J-1)';                      % (J-1)*2
c_mu = c(:,1);
c_se = c(:,2);

seed = 1;
prob_simu = zeros(I,1);
Respondent = Respondent_mat(:,1);
threshold = 4;

for i = 1:I
    rng(seed);
%     trun_randn = randn(1,(K+J-1)*D);            % w: Still? User objective function returned NaN; trying a new point
%     if sum(abs(trun_randn)>threshold)>=1
%         trun_randn(abs(trun_randn)>threshold) = [];
%     end
%     while length(trun_randn)<(K+J-1)*D
%         seed = seed +1;
%         trun_randn = [trun_randn randn(1,(K+J-1)*D-length(trun_randn))];
%         if sum(abs(trun_randn)>threshold)>=1
%             trun_randn(abs(trun_randn)>threshold) = [];
%         end
%     end
%     trun_randn = reshape(trun_randn,K+J-1,D);
%     idraws = repmat(mu,1,D) + trun_randn(1:K,:).*repmat(exp(se),1,D);         % K*D        nVar*nDraws
%     cdraws = repmat(c_mu,1,D)+trun_randn(K+1:K+J-1,:).*repmat(exp(c_se),1,D);     %(J-1)*D
    constrain_randn = randn(1,(K+J-1)*D*2);          % over sample by 1 time
    constrain_randn(abs(constrain_randn)>threshold) = [];
    constrain_randn = constrain_randn(1,1:(K+J-1)*D);
    constrain_randn = reshape(constrain_randn,K+J-1,D);
    idraws = repmat(mu,1,D) + constrain_randn(1:K,:).*repmat(exp(se),1,D);         % K*D        nVar*nDraws
    cdraws = repmat(c_mu,1,D)+constrain_randn(K+1:K+J-1,:).*repmat(exp(c_se),1,D);     %(J-1)*D
%     idraws = repmat(mu,1,D) + randn(K,D).*repmat(exp(se),1,D);        % K*D        nVar*nDraws
%     cdraws = repmat(c_mu,1,D) + randn(J-1,D).*repmat(exp(c_se),1,D);     %(J-1)*D
    const = [cdraws;zeros(1,D)];                                 % J*D     last alternative as base line       
    [r c p] = find(Respondent==i);
    iIV = IV(r,:,:);                               % T*J*K
    iTime = length(r);
    ichoices = choice_dv(r,:);             % T*J
    % p_seq = ones(1,D);
    sum_pt = zeros(1,D);
    for t = 1:iTime
        expU_t = exp(squeeze(iIV(t,:,:))*idraws + const);            % J*D
        % expU_t = exp(util(X, bs))
        p_chosen_t = expU_t(ichoices(t,:)==1,:)./sum(expU_t,1);
        sum_pt = sum_pt + log(p_chosen_t);
        % p_seq = p_seq.*p_chosen_t;
    end
    p_seq = exp(sum_pt);
    prob_simu(i,1) = mean(p_seq);
    seed = seed+1;              % w: ????????????
end

LL = -sum(log(prob_simu));

end