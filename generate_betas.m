function [beta1, beta2, beta3, beta4, beta6] = generate_betas(h, I, J)

beta1_user = exp(randn(I,1));
beta2_user = exp(randn(I,1));
beta3_user = exp(randn(I,1));
beta4_user = exp(randn(I,1));
beta6_user = exp(randn(I,1));

beta1_alt = exp(randn(1,J));
beta2_alt = exp(randn(1,J));
beta3_alt = exp(randn(1,J));
beta4_alt = exp(randn(1,J));
beta6_alt = exp(randn(1,J));

beta1 = beta1_user*beta1_alt.*h(1);
beta2 = beta2_user*beta2_alt.*h(2);
beta3 = beta3_user*beta3_alt.*h(3);
beta4 = beta4_user*beta4_alt.*h(4);
beta6 = beta6_user*beta6_alt.*h(6);

end