% function [LL, gr, H] = band_bi_ll_i2(b,IVs,choice_dv, innov_X, explor_X, week_IV, innov_WD_multip, explor_WD_multip)
function [LL, g, H] = band_bi_ll_trend_reverse_N_age_diverse_ms(b,IVs,trend_hat, band_age, topics_count, X_N, None0s_X_N, S, choice_dv)
global dummies se

I = size(choice_dv,1);
% K = size(IVs,2);

const = b(1);

% FV = IVs*bs;
b_basic = b(4:13)';

week_IV_temp = 100*gampdf(IVs(:,2),exp(b(2)),exp(b(3)));          
week_IV=week_IV_temp.*IVs(:,3);

val_pdf = 100*normpdf(None0s_X_N(:,3),0,b(14));
X_N = sparse(None0s_X_N(:,1),None0s_X_N(:,2),val_pdf,50991509,6222);
IV_N_S = X_N*S.^exp(b(15));

FV = [IVs(:,1) trend_hat week_IV  band_age  topics_count...
    band_age.*week_IV  topics_count.*week_IV   band_age.*topics_count...
    band_age.*topics_count.*week_IV    IV_N_S]*b_basic;          % with both trend and BP as controls

% exp_util = exp(const+FV);          % utility of choosing the product
% prob=exp_util./(1+exp_util);
exp_util = exp(-(const+FV));         % this is now the utility of the external good
prob=1./(1+exp_util);                % this is still the probability of choosing the product
pmat = [prob 1-prob]; 
pmat = pmat.*choice_dv;
[r c p] = find(pmat);                                             % I*1
LL = -sum(log(p));                                                % 1*1

% g(1) = sum(prob-choice_dv(:,1));
% g(4:13) = [IVs(:,1)  trend_hat  week_IV  band_age  topics_count...
%     band_age.*week_IV  topics_count.*week_IV   band_age.*topics_count...
%     band_age.*topics_count.*week_IV   IV_N_S]'*(prob-choice_dv(:,1));
% gamma_derivative_alpha = gamma_derivative_wrt_alpha(IVs(:,2), b(2), b(3));
% g(2) = 100*b(6).*IVs(:,3)'.*gamma_derivative_alpha'*(prob-choice_dv(:,1))+...
%     100*b(9).*band_age'.*IVs(:,3)'.*gamma_derivative_alpha'*(prob-choice_dv(:,1))+...
%     100*b(10).*topics_count'.*IVs(:,3)'.*gamma_derivative_alpha'*(prob-choice_dv(:,1))+...
%     100*b(12).*band_age'.*topics_count'.*IVs(:,3)'.*gamma_derivative_alpha'*(prob-choice_dv(:,1));
% gamma_derivative_beta = gamma_derivative_wrt_beta(IVs(:,2), b(2), b(3));
% g(3) = 100*b(6).*IVs(:,3)'.*gamma_derivative_beta'*(prob-choice_dv(:,1))+...
%     100*b(9).*band_age'.*IVs(:,3)'.*gamma_derivative_beta'*(prob-choice_dv(:,1))+...
%     100*b(10).*topics_count'.*IVs(:,3)'.*gamma_derivative_beta'*(prob-choice_dv(:,1))+...
%     100*b(12).*band_age'.*topics_count'.*IVs(:,3)'.*gamma_derivative_beta'*(prob-choice_dv(:,1));
% normal_derivative_wrt_sigmas = normal_derivative(None0s_X_N(:,3), 0, b(14)); 
% X_N_derivative = sparse(None0s_X_N(:,1),None0s_X_N(:,2),100.*normal_derivative_wrt_sigmas,50991509,6222);
% g(14) = b(13).*(X_N_derivative*S.^exp(b(15)))'*(prob-choice_dv(:,1));
% g(15) = b(13).*(X_N*similarity_derivative_wrt_coef(S, b(15)))'*(prob-choice_dv(:,1));

end