function g = gamma_derivative_wrt_beta(x, alpha, beta)           % normal_derivative wrt sigmas

% syms f alpha beta x
% f = (x^(alpha - 1)*exp(-beta)^alpha*exp(-x*exp(-beta)))/gamma(alpha)
% g = diff(f,beta)

g = (x.*x.^(alpha - 1).*exp(-beta).^alpha.*exp(-beta).*exp(-x.*exp(-beta)))./gamma(alpha)...
    - (alpha.*x.^(alpha - 1).*exp(-beta).^(alpha - 1).*exp(-beta).*exp(-x.*exp(-beta)))./gamma(alpha);

end