clc
clear
global I IV choice_dv J K Respondent_mat

rng(0);

J = 5
I = 200
T = 50
K = 0
% unobs0 = ones(1,J).*2;
% unobs1 = ones(1,J);

unobs0 = [1 2 3 4 5]
unobs1 = [0.5 1.5 1 2 1]
% unobs0 = [1 2 3 4]
% unobs1 = [0.5 1.5 1 2]

% unobs0 = 3
% unobs1 = 2
const_alt = zeros(I*T,J);
for j = 1:J                  %  CHANGED !
    mu_j = unobs0(j) + randn(I,1).*unobs1(j);
    const_alt(:,j) = reshape(repmat(mu_j,1,T)',I*T,1);
end

features = randn(I*T,J,K);
features=features-mean(mean(mean(features)));
IV = features;

% beta0 = ones(1,K);                       % w: E: beta is cost or price
% beta1 = ones(1,K);
beta0 = [3 7]
beta1 = [1 1.5]
beta_alt = zeros(I*T,K);                 % IT*K
for k = 1:K
    muj = beta0(k) + randn(I,1).*beta1(k);                 % I*1
    beta_alt(:,k) = reshape(repmat(muj,1,T)',I*T,1);        % IT*1
end

rep_beta_alt = permute(repmat(beta_alt,[1,1,J]),[1 3 2]);
% exp_utility = exp(const_alt + sum(features.*rep_beta_alt,3));          % IT*J      w!!!!!
exp_utility = exp(const_alt);
prob = exp_utility./repmat(sum(exp_utility,2),1,J);               % IT*J
prob=cumsum(prob')';
draw_for_choice=rand(I*T,1);
draw_for_choice=repmat(draw_for_choice,1,J);
choice=prob<draw_for_choice;
choice=sum(choice,2)+1;

choicemat = repmat(choice,1,J);
testmat = repmat(1:J,I*T,1);
choice_dv = choicemat==testmat;

individuals = 1:I;
Respondents = reshape(repmat(individuals',1,T)',I*T,1);
Respondent_mat = repmat(Respondents,1,J);

tic

b0 = zeros(1,2*(J-1+K))

% options = optimset('LargeScale','off','GradObj','off','Hessian','off','display','iter')
options = optimset('LargeScale','off','GradObj','off','Hessian','off','display','iter', 'MaxIter',1e4, 'MaxFunEvals', 1e5)         % LargeScale off is quasi-Newton method in optimset

[beta_0, fval,exitflag,output,grad,hessian] = fminunc(@simp_ll_test4,b0,options);

se = sqrt(diag(inv(hessian)))

betas = reshape(beta_0,2,[])';
betas(:,2) = exp(betas(:,2))

toc