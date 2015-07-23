function LL = mnl_ll(beta_0)
global I IV choice_dv J K Respondent_mat cswitch

D = 100;                                        % number of draws

if cswitch == 1
    % bb = beta_0(2*J-1:2*(J-1+K));
    bb = beta_0(2*J+1:2*(J+K));
    b = reshape(bb,2,K)';                            % K*2
    mu = b(:,1);
    se = b(:,2);
    
    % bc = beta_0(1:2*J-2);
    % c = reshape(bc,2,J-1)';                      % (J-1)*2
    bc = beta_0(1:2*J);
    c = reshape(bc,2,J)';                      % (J-1)*2
    c_mu = c(:,1);
    c_se = c(:,2);
    c_mu(J) = 0;
    
    % c_mu = zeros(size(c_mu));                                 % to fix all constants to 0
    % c_se = zeros(size(c_se));
    % % c_mu(1) = 1;                            % fixed
else
    bb = beta_0;
    b = reshape(bb,2,K)';                            % K*2
    mu = b(:,1);
    se = b(:,2);
end

seed = 14239;
% seed = 100;
prob_simu = zeros(I,1);
Respondent = Respondent_mat(:,1);
threshold = 4;

% draws_collect = zeros(I,(K+J)*D);
for i = 1:I
    rng(seed);
%     constrain_randn = randn(1,(K+J-1)*D*2);          % over sample by 1 time
%     constrain_randn(abs(constrain_randn)>threshold) = [];
%     constrain_randn = constrain_randn(1,1:(K+J-1)*D);
%     constrain_randn = reshape(constrain_randn,K+J-1,D);                           % (K+J-1)*D
%     idraws = repmat(mu,1,D) + constrain_randn(1:K,:).*repmat(exp(se),1,D);         % K*D        nVar*nDraws
%     cdraws = repmat(c_mu,1,D)+constrain_randn(K+1:K+J-1,:).*repmat(exp(c_se),1,D);     %(J-1)*D
    constrain_randn = randn(1,(K+J)*D*1.5);          % over sample by 0.5 time
    constrain_randn(abs(constrain_randn)>threshold) = [];
    constrain_randn = constrain_randn(1,1:(K+J)*D);
%    draws_collect(i,:) = constrain_randn;
    constrain_randn = reshape(constrain_randn,K+J,D);                           % (K+J)*D
    idraws = repmat(mu,1,D) + constrain_randn(1:K,:).*repmat(exp(se),1,D);         % K*D        nVar*nDraws
    if cswitch == 1
        cdraws = repmat(c_mu,1,D)+constrain_randn(K+1:K+J,:).*repmat(exp(c_se),1,D);     %(J-1)*D
    else
        cdraws = 0;
    end
    %     idraws = repmat(mu,1,D) + randn(K,D).*repmat(exp(se),1,D);        % K*D        nVar*nDraws
%     cdraws = repmat(c_mu,1,D) + randn(J-1,D).*repmat(exp(c_se),1,D);     %(J-1)*D
    % const = [cdraws;zeros(1,D)];                                 % J*D     last alternative as base line 
    const = cdraws;
    [r c p] = find(Respondent==i);
    iIV = IV(r,:,:);                               % T*J*K
    iTime = length(r);
    ichoices = choice_dv(r,:);             % T*J
    % p_seq = ones(1,D);
    sum_pt = zeros(1,D);
    for t = 1:iTime
        if K == 1
           expU_t = exp(iIV(t,:,:)'*idraws + const);
        else
           expU_t = exp(squeeze(iIV(t,:,:))*idraws + const);            % J*D
        end
        % expU_t = exp(util(X, bs))
        p_chosen_t = expU_t(ichoices(t,:)==1,:)./sum(expU_t,1);      % 1*D
        sum_pt = sum_pt + log(p_chosen_t);
        % p_seq = p_seq.*p_chosen_t;
    end
    p_seq = exp(sum_pt);
    prob_simu(i,1) = mean(p_seq);
    seed = seed+1;              % w: ????????????
end

LL = -sum(log(prob_simu));

end