function f = loglike(vmus, l1, l2, hedonic) 
beta0=vmus(1);
beta1=vmus(2);
beta2=vmus(3);
% PROB = 1./(1+exp(-beta0 - beta1.*l1 - beta2.*l2));
% aux=hedonic.*PROB + (1-hedonic).*(1-PROB);
PROB = exp(beta0 + beta1.*l1 + beta2.*l2)./(1+exp(beta0 + beta1.*l1 + beta2.*l2));
aux=hedonic.*PROB + (1-hedonic).*(1-PROB);
f=-sum(log(aux)); 