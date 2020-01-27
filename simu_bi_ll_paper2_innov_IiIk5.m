% function [LL, gr, H] = band_bi_ll_i2(b,IVs,choice_dv, innov_X, explor_X, week_IV, innov_WD_multip, explor_WD_multip)
function [LL, gr, H] = simu_bi_ll_paper2_innov_IiIk5(b,IVs,trend_hat, pop, dummy_agg_SI, mat_sparse_8088, m_innov, f_innov, Im_17617085, Ik_17617085, cosine_similarity_scores, None0s_X_N, S, choice_dv,rm_choice_Innovf)
global dummies se

I = size(choice_dv,1);
% K = size(IVs,2);

const = b(1);

% FV = IVs*bs;
b_basic = b(2:8)';

% week_IV = dummy_agg_SI;
% week_IV_innov = dummy_agg_SI_innov;

Aijkt = mat_sparse_8088;
Sik = cosine_similarity_scores;

val_pdf = 100*normpdf(None0s_X_N(:,3),0,b(9));
X_N = sparse(None0s_X_N(:,1),None0s_X_N(:,2),val_pdf,17617085,6222);
IV_N_S = X_N*S.^exp(b(10));

pop = pop./1000;

f_innov = f_innov*10;
m_innov = m_innov*10;
Aijkt = Aijkt*10;

% Sik_power = Sik.^exp(b(11));
% Sik_power = Sik;
Sik_power = ones(8088,1);
% IV_N_S = IV_N_S*10;

% Im_17617085 = Im_17617085./10;
% Ik_17617085 = Ik_17617085./10;

% FV = [IVs(:,1)  trend_hat  week_IV  band_age  topics_count...
%     band_age.*week_IV  topics_count.*week_IV   band_age.*topics_count...
%     band_age.*topics_count.*week_IV    IV_N_S]*b_basic;          % with both trend and BP as controls
% FV = [trend_hat  week_IV  week_IV.^2   pop   week_IV.*pop      week_IV.^2.*pop     week_IV_innov    IV_N_S]*b_basic;          % with both trend and BP as controls
FV = [Im_17617085    Ik_17617085     Aijkt*Sik_power    Aijkt*(m_innov.*Sik_power)    Aijkt*(f_innov.*Sik_power.*rm_choice_Innovf)    Aijkt*(f_innov.*m_innov.*Sik_power.*rm_choice_Innovf)     IV_N_S]*b_basic;

% week_IV.^2*pop

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