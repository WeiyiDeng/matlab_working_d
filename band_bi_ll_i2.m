% function [LL, gr, H] = band_bi_ll_i2(b,IVs,choice_dv, innov_X, explor_X, week_IV, innov_WD_multip, explor_WD_multip)
function [LL, gr, H] = band_bi_ll_i2(b,IVs,choice_dv, innov_X, explor_X)
global dummies se

% b = (diag([1e-1,1e-1,1e2])*b')';    % w: this does not yield the right result

% const = repmat(b(1),I*T,1);

I = size(choice_dv,1);
% K = size(IVs,2);

const = b(1);

bs = b(4:end)';
% bs = b(2:end)';
% bs = b(2:end)';
% bs = b(3:end)';
% bs = b(2:1+size(IVs,2))';
% FV = IVs*bs;

% if k<2, when x = 0, fk(x)= inf
% documentation of chi2pdf requires the degree of freedom parameter k must be positive integers
% chi-square dist can be considered as a generalization of gamma distribution, which can be extended to non integer values
% week_IV = IVs(:,3).*chi2pdf(IVs(:,2),b(2));     
% week_IV = IVs(:,3).*gampdf(IVs(:,2),exp(b(2)),exp(b(3)));         % exp() to change estimated beta to larger than zero

% % week_IV = IVs(:,3).*gampdf(IVs(:,2),b(2),b(3));          % no need to multiply IVs(:,3) since will set f(x<=0) = 0 
% week_IV = gampdf(IVs(:,2),b(2),b(3));   
% % week_IV = gampdf(IVs(:,2),1.1,20);    % in the case when the first parameter of gamma is larger than 1, at x=0 the prob will be 0 anyway           
% % week_IV = IVs(:,3).*exppdf(IVs(:,2),b(2));
% week_IV(IVs(:,2)<1)=0;                  % when x equal to (or smaller than) 0 the variable is set to zero, not to be estimated 

% triang_distr = @(x) (b(2)-x)*2/((b(2)-1)*(b(2)-1));
% week_IV = triang_distr(IVs(:,2));
% week_IV(IVs(:,2)<1 | IVs(:,2)>b(2))=0;

% WD_IV = IVs(:,2);
% gamma_trans = zeros(length(WD_IV),1);
% parfor i = 1:length(WD_IV)
%     gamma_trans(i) = gampdf(WD_IV(i),b(2),b(3));            % w: cpu stucked here, seems to be slower with separate indices. Why?
% end
% week_IV = IVs(:,3).*gamma_trans;

% week_IV = IVs(:,3).*gampdf(IVs(:,2),1,2);
% FV = [IVs(:,1) week_IV]*[bs; 0.07];

week_IV = gampdf(IVs(:,2),b(2),b(3));          
week_IV(IVs(:,2)<1)=0;      

innov_WD_multip = zeros(size(innov_X));
for i = 1:size(innov_X,2);
    innov_WD_multip(:,i) = innov_X(:,i).*week_IV;
end

% explor_WD_multip = zeros(size(explor_X));
% for j = 1:size(explor_X,2);
%     explor_WD_multip(:,j) = explor_X(:,j).*week_IV;
% end
explor_WD_multip = [];

% FV = [IVs(:,1) week_IV]*[3.0426   12.7745]'; 
FV = [IVs(:,1) week_IV innov_X explor_X innov_WD_multip explor_WD_multip]*bs;
% new_FV = [innov_X explor_X innov_WD_multip explor_WD_multip]*bs;
% FV = FV + new_FV;
% clearvars week_IV WD_IV gamma_trans

% bs_d = b(2+size(IVs,2):end);
% D = zeros(I,1);
% for i = 1:size(dummies,1)
%     D(i) = sum(dummies(i,:).*bs_d);       % avoid matrix multiplication in sparse matrix
% end
%     
% FV = FV+D;

% exp_util = exp(const+FV);          % utility of choosing the product
% prob=exp_util./(1+exp_util);
exp_util = exp(-(const+FV));         % this is now the utility of the external good
prob=1./(1+exp_util);                % this is still the probability of choosing the product
pmat = [prob 1-prob]; 
pmat = pmat.*choice_dv;
[r c p] = find(pmat);                                             % I*1
LL = -sum(log(p));                                                % 1*1

% format long g
% LL

% d = choice_dv(:,1)-prob;                                    % I*1                     
% gr = d'*[ones(I,1) IVs];                          % 1*k
% gr = -gr;                                                   % w: why reverse ??
% 
% H = (repmat(prob.*(1-prob),1,size(IVs,2)+1).*[ones(I,1) IVs])'*[ones(I,1) IVs];

% d = prob - choice_dv(:,1);                                % I*1        w: reverse?    
% Gt = repmat(d,1,size(b,2)).*[ones(size(IVs,1),1) IVs];    % I*k
% gr = sum(Gt,1);

% small_value = (2.*rand(1,3)-1)*0.0000005;
% gr = sum(Gt,1).*(1+small_value);

% d = pmat-choice_dv;                                % I*2
% Gt = zeros(size(IVs,1), size(IVs,2)+1);
% Gt(:,1) = sum([zeros(size(IVs,1),1) -ones(size(IVs,1),1)].*d, 2);
% for i = 1:size(IVs,2)
%     Gt(:,i+1) = sum([zeros(size(IVs,1),1) -IVs(:,i)].*d, 2);
% end
% gr = sum(Gt,1);

% % w: hessian not working, why ??
% H = inv(Gt'*Gt);                          % H is actually the inverse of B here (B is approaching -Ht asymptotically)
% disp(gr);
% se = sqrt(diag(H));                      % H is the covariance matrix here ?

%  p_all = utility_all./repmat(sum(utility_all,2),1,J);
%  d = p_all-choice_dv;                     % beware the sequence of what minus what !!
%  
%  Gt = zeros(I,(2*(J-1)+1));                         % number of variables
% %  Gt(:,1) = sum(ones(I,1).*d(:,1), 2);
% %  Gt(:,2) = sum(ones(I,1).*d(:,2), 2);
% %  Gt(:,3) = sum(ones(I,1).*d(:,3), 2);
% %  Gt(:,4) = sum(ones(I,1).*d(:,4), 2);
% for j = 1:(J-1)
%     Gt(:,j) = sum(ones(I,1).*d(:,j), 2);
%     Gt(:,j+1+J-1) = sum(demograph.*d(:,j), 2);
% end
% Gt(:,5) = sum(IV2.*d, 2);
% %  Gt(:,6) = sum(demograph.*d(:,1), 2);
% %  Gt(:,7) = sum(demograph.*d(:,2), 2);
% %  Gt(:,8) = sum(demograph.*d(:,3), 2);
% %  Gt(:,9) = sum(demograph.*d(:,4), 2);
%  
%  gr = sum(Gt,1);
%  
%  H = inv(Gt'*Gt);                          % H is actually the inverse of B here (B is approaching -Ht asymptotically) 
%  disp(gr);
%  se = sqrt(diag(H));                      % H is the covariance matrix here ? 
%  % disp(se);
%  % fprintf('standard error %d ',se);

end