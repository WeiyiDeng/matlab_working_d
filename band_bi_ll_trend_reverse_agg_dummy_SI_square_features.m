% function [LL, gr, H] = band_bi_ll_i2(b,IVs,choice_dv, innov_X, explor_X, week_IV, innov_WD_multip, explor_WD_multip)
function [LL, gr, H] = band_bi_ll_trend_reverse_agg_dummy_SI_square_features(b,IVs,trend_hat, pop, topics_count, band_age, tracks_count, dummy_agg_SI, None0s_X_N, S, choice_dv)
global dummies se

I = size(choice_dv,1);
% K = size(IVs,2);

const = b(1);

% FV = IVs*bs;
b_basic = b(2:18)';

week_IV = dummy_agg_SI;

val_pdf = 100*normpdf(None0s_X_N(:,3),0,b(19));
X_N = sparse(None0s_X_N(:,1),None0s_X_N(:,2),val_pdf,17617085,6222);
IV_N_S = X_N*S.^exp(b(20));

pop = pop./100;
topics_count = topics_count./100;
band_age = band_age./1000;
tracks_count = tracks_count./1000;

% FV = [IVs(:,1)  trend_hat  week_IV  band_age  topics_count...
%     band_age.*week_IV  topics_count.*week_IV   band_age.*topics_count...
%     band_age.*topics_count.*week_IV    IV_N_S]*b_basic;          % with both trend and BP as controls
FV = [IVs(:,1)  trend_hat  week_IV  week_IV.^2   pop   topics_count     band_age     tracks_count...
    week_IV.*pop      week_IV.^2.*pop     week_IV.*topics_count      week_IV.^2.*topics_count...
    week_IV.*band_age      week_IV.^2.*band_age     week_IV.*tracks_count      week_IV.^2.*tracks_count   IV_N_S]*b_basic;          % with both trend and BP as controls

% week_IV.^2.*pop

% exp_util = exp(const+FV);          % utility of choosing the product
% prob=exp_util./(1+exp_util);
exp_util = exp(-(const+FV));         % this is now the utility of the external good
prob=1./(1+exp_util);                % this is still the probability of choosing the product
pmat = [prob 1-prob]; 
pmat = pmat.*choice_dv;
[r c p] = find(pmat);                                             % I*1
LL = -sum(log(p));                                                % 1*1

% assert(det(H)~=0, 'Hessian is singular');

end