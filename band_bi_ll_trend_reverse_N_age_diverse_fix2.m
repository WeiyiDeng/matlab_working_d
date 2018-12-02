% function [LL, gr, H] = band_bi_ll_i2(b,IVs,choice_dv, innov_X, explor_X, week_IV, innov_WD_multip, explor_WD_multip)
function [LL, gr, H] = band_bi_ll_trend_reverse_N_age_diverse_fix2(b,IVs,trend_hat, band_age, topics_count, X_N, None0s_X_N, S, choice_dv, beta_fix)
global dummies se

I = size(choice_dv,1);
% K = size(IVs,2);

const = beta_fix(1);

% FV = IVs*bs;
b_basic = beta_fix(4:9)';

week_IV = 100*gampdf(IVs(:,2),exp(beta_fix(2)),exp(beta_fix(3)));          
week_IV(IVs(:,3)==0)=0;

val_pdf = 100*normpdf(None0s_X_N(:,3),0,exp(beta_fix(10)));
X_N = sparse(None0s_X_N(:,1),None0s_X_N(:,2),val_pdf,50991509,6222);
IV_N_S = X_N*S.^exp(beta_fix(11));

FV = [IVs(:,1) trend_hat week_IV  band_age  topics_count  IV_N_S...
    band_age.*week_IV  topics_count.*week_IV   band_age.*topics_count...
    band_age.*topics_count.*week_IV]*[b_basic; b'];          % with both trend and BP as controls

% exp_util = exp(const+FV);          % utility of choosing the product
% prob=exp_util./(1+exp_util);
exp_util = exp(-(const+FV));         % this is now the utility of the external good
prob=1./(1+exp_util);                % this is still the probability of choosing the product
pmat = [prob 1-prob]; 
pmat = pmat.*choice_dv;
[r c p] = find(pmat);                                             % I*1
LL = -sum(log(p));                                                % 1*1

end