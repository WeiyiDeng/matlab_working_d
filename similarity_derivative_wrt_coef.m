function g = similarity_derivative_wrt_coef(S, alpha)           % normal_derivative wrt alpha

% syms f alpha S
% f = S^exp(alpha)
% g = diff(f,alpha)

g = S.^exp(alpha).*exp(alpha).*log(S);

end