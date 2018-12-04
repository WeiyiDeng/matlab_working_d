function g = gamma_derivative_wrt_alpha(x, alpha, beta)           % normal_derivative wrt alpha
%w% gamma_derivative_wrt_alpha not defined for x=0(returns NAN)

% syms f alpha beta x
% f = ((1/exp(beta))^exp(alpha)*x^(exp(alpha)-1)*exp(-x/exp(beta)))/gamma(exp(alpha))
% g = diff(f,alpha)

g = (x.^(exp(alpha) - 1).*exp(-beta).^exp(alpha).*exp(-x.*exp(-beta)).*exp(alpha).*log(x))./gamma(exp(alpha))...
    - (x.^(exp(alpha) - 1).*exp(-beta).^exp(alpha).*psi(exp(alpha)).*exp(-x.*exp(-beta)).*exp(alpha))./gamma(exp(alpha))...
    + (x.^(exp(alpha) - 1).*exp(-beta).^exp(alpha).*exp(-x.*exp(-beta)).*exp(alpha).*log(exp(-beta)))./gamma(exp(alpha));

end