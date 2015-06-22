function [iv, dv, respond] = simu(J,I,T,K)

rng(0);

unobs0 = (1:J).*2;           % means
unobs1 = ones(1,J);          % variances

const_alt = zeros(I*T,J);
for j = 1:J
    mu_j = unobs0(j) + randn(I,1).*unobs1(j);
    const_alt(:,j) = reshape(repmat(mu_j,1,T)',I*T,1);
end

features = randn(I*T,J,K);
features=features-mean(mean(mean(features)));
iv = features;

beta0 = ones(1,K).*2;
beta1 = ones(1,K);

beta_alt = zeros(I*T,K);                 % IT*K
for k = 1:K
    muj = beta0(k) + randn(I,1).*beta1(k);                 % I*1
    beta_alt(:,k) = reshape(repmat(muj,1,T)',I*T,1);        % IT*1
end

rep_beta_alt = permute(repmat(beta_alt,[1,1,J]),[1 3 2]);
exp_utility = exp(const_alt + sum(features.*rep_beta_alt,3));          % IT*J      w!!!!!
% exp_utility = exp(util(X, bs))
prob = exp_utility./repmat(sum(exp_utility,2),1,J);               % IT*J
prob=cumsum(prob')';
draw_for_choice=rand(I*T,1);
draw_for_choice=repmat(draw_for_choice,1,J);
choice=prob<draw_for_choice;
choice=sum(choice,2)+1;

choicemat = repmat(choice,1,J);
testmat = repmat(1:J,I*T,1);
dv = choicemat==testmat;

individuals = 1:I;
Respondents = reshape(repmat(individuals',1,T)',I*T,1);
respond = repmat(Respondents,1,J);