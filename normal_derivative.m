function g = normal_derivative(x, mu, sigmas)           % normal_derivative wrt sigmas

% syms f mu sigmas pai x
% f=1/(sqrt(2.*pai*exp(sigmas))).*exp(-(x-mu)^2./(2.*exp(sigmas)));
% g = diff(f,sigmas)

g = (2.^(1/2).*exp(-(exp(-sigmas).*(mu - x).^2)./2).*exp(-sigmas).*(mu - x).^2)./(4.*(pi.*exp(sigmas)).^(1/2))...
    - (2.^(1/2).*pi.*exp(-(exp(-sigmas).*(mu - x).^2)./2).*exp(sigmas))./(4*(pi.*exp(sigmas)).^(3/2));

end