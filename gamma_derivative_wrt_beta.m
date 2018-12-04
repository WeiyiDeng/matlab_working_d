function g = gamma_derivative_wrt_beta(x, alpha, beta)           % normal_derivative wrt beta

% % syms f alpha beta x
% % f = ((1/exp(beta))^alpha*x^(alpha-1)*exp(-x/exp(beta)))/gamma(alpha)
% % g = diff(f,beta)
% 
% g = (x.*x.^(alpha - 1).*exp(-beta).^alpha.*exp(-beta).*exp(-x.*exp(-beta)))./gamma(alpha)...
%     - (alpha.*x.^(alpha - 1).*exp(-beta).^(alpha - 1).*exp(-beta).*exp(-x.*exp(-beta)))./gamma(alpha);

% syms f alpha beta x
% f = ((1/exp(beta))^exp(alpha)*x^(exp(alpha)-1)*exp(-x/exp(beta)))/gamma(exp(alpha))
% g = diff(f,beta)

g = (x.*x.^(exp(alpha) - 1).*exp(-beta).^exp(alpha).*exp(-beta).*exp(-x.*exp(-beta)))./gamma(exp(alpha))...
    - (x.^(exp(alpha) - 1).*exp(-beta).^(exp(alpha) - 1).*exp(-beta).*exp(-x.*exp(-beta)).*exp(alpha))./gamma(exp(alpha));

end