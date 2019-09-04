bar = length(f_innov_percent);
prob = zeros(bar,bar);
for m = 1:bar
    ind_m = m;
    for k = 1:bar
        ind_k = k;
        % FV = [Im_percent(ind)    Ik_percent(ind)     SI_percent(ind)    SI_percent(ind)*m_innov_percent(ind)    SI_percent(ind)*f_innov_percent(ind)    SI_percent(ind)*m_innov_percent(ind)*f_innov_percent(ind)     N_percent(ind)]*b_basic;
        FV = [Im_percent(ind_m)    Ik_percent(ind_k)     SI_percent    SI_percent*m_innov_percent(ind_m)    SI_percent*f_innov_percent(ind_k)    SI_percent*m_innov_percent(ind_m)*f_innov_percent(ind_k)     N_percent]*b_basic;
        exp_util = exp(-(const+FV));         % this is now the utility of the external good
        prob(m,k)=1./(1+exp_util)                % this is still the probability of choosing the product
    end
end
surf(prob)