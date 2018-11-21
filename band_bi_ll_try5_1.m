% function [LL, gr, H] = band_bi_ll_i2(b,IVs,choice_dv, innov_X, explor_X, week_IV, innov_WD_multip, explor_WD_multip)
function [LL, gr, H] = band_bi_ll_try5_1(b,IVs,choice_dv)
global dummies se

const = b(1);

FV_contrl = IVs(:,1:end-4)*b(2:end-2)';
FV_X = IVs(:,end-3)+IVs(:,end-2).*(1/exp(b(end))^1)+IVs(:,end-1).*(1/exp(b(end))^2)+IVs(:,end).*(1/exp(b(end))^3);

exp_util = exp(-(const+FV_contrl+b(end-1)*FV_X));         % this is now the utility of the external good
prob=1./(1+exp_util);                % this is still the probability of choosing the product
pmat = [prob 1-prob]; 
temp = pmat.*choice_dv;
[r c p] = find(temp);                                             % I*1
LL = -sum(log(p));                                                % 1*1

end