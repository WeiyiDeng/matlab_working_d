function [betas, se, real_bs] = run_mnl(data, Xlocation, ylocation, Ilocation, b0)
global I IV choice_dv J N K Respondent_mat PARAMTOL

if any(data) == 1
    N = size(data,1);
    J = length(ylocation);
    I = data(N,Ilocation(1));
    K = length(Xlocation)/J;
    IV = zeros(N,J,K);
    data_iv = data(:,Xlocation);
    for k = 1:K
        IV(:,:,k) = data_iv(:,((k-1)*J+1):k*J);
    end
    Respondent_mat = data(:,Ilocation);
    choice_dv = data(:,ylocation);
    real_bs = [];
else
    if isempty(ylocation) == 1
        [N, J, K] = size(data);
        T = Xlocation;
        I = N/T;
        [iv, dv, respond, real_bs] = simu(J,I,T,K);
        IV = iv;
        Respondent_mat = respond;
        choice_dv = dv;
    else
        [N, J, K] = size(data);
        T = Xlocation;
        I = N/T;
        [iv, dv, respond, real_bs] = debug(J,I,T,K);
        IV = iv;
        Respondent_mat = respond;
        choice_dv = dv;
    end
end

options = optimset('LargeScale','off','GradObj','off','Hessian','off','display','iter', 'MaxIter',1e4, 'TolX',PARAMTOL,'MaxFunEvals', 1e5)

[beta_0, fval,exitflag,output,grad,hessian] = fminunc(@mnl_ll,b0,options);

se = sqrt(diag(inv(hessian)))

betas = reshape(beta_0,2,[])';
betas(:,2) = exp(betas(:,2)) 

end