clc
% clear
clear all

% diary('wwdiary.txt')

% date_=clock;
% resultsfilename=['Results/r_Results_' num2str(date_(1)) '_' num2str(date_(2)) '_' num2str(date_(3)) '-' num2str(date_(4))  '_' num2str(date_(5)) '.txt'];
% diary(resultsfilename);

% load('matp_b.mat');
% load('matp.mat');
% load('innov_contin.mat');
% load('explor_contin.mat');
load('matpstd.mat');

% % w: draw a small random sample from the obs 
% temp1=rand(size(matp,1),1)>.985;     
% matp=matp(  matp(:,5)==1 | (matp(:,5)==0 &  temp1 ), :);
% matp=matp(rand(size(matp,1),1)>.95, :);
% sum(matp(:,5)==1)

% num = 200
% rand_vector = rand(num,1);
% y_simu = zeros(num,1);
% y_simu(rand_vector>0.5) = 1;
% 
% matp = zeros(num,7);
% matp(:,5) = y_simu;

% mat1 = csvread('imat4474.csv');
% load('member_dummies.mat');
% load('member_dummies_week_d.mat');

% matp(:,6) = matp(:,6)-mean(matp(:,6));
% matp(:,7) = matp(:,7)-mean(matp(:,7));
% matp(:,6) = (matp(:,6)-mean(matp(:,6)))./std(matp(:,6));
% matp(:,7) = (matp(:,7)-mean(matp(:,7)))./std(matp(:,7));

%%
% load('matpstd.mat');
% 
% X = matp(:,[6 7]);
% y = matp(:,5);
% 
% IVs = X;
% choice_dv = [y 1-y];
% 
% week_IV = 100*gampdf(IVs(:,2),0.8375,27.4694);              % w: NOTICE new lines here for fixed gamma parameters
% week_IV(IVs(:,2)<1)=0; 
% 
% prev_FV_basic = [IVs(:,1) week_IV]*[0.3870    0.1020]';
% 
% clearvars matp IVs
% 
% load('innov_contin_std.mat');
% % load('explor_contin_std.mat');
% innov_X = innov_contin(:,1:4);
% 
% clearvars innov_contin
% 
% innov_WD_multip = zeros(size(innov_X));
% for i = 1:size(innov_X,2);
%     innov_WD_multip(:,i) = innov_X(:,i).*week_IV;
% end
% % innov_WD_multip = [];
% 
% prev_FV_innov = [innov_X innov_WD_multip]*[-0.0509    0.0281   -0.2486    0.2183  0.0102    0.0140    0.0414   -0.0358]';
% 
% clearvars innov_X innov_WD_multip
% 
% load('explor_contin_std.mat');
% explor_X = explor_contin(:,1:4);
% 
% clearvars explor_contin
% 
% explor_WD_multip = zeros(size(explor_X));
% for j = 1:size(explor_X,2);
%     explor_WD_multip(:,j) = explor_X(:,j).*week_IV;
% end
% % explor_WD_multip = [];
% 
% prev_FV_explor = [explor_X explor_WD_multip]*[0.1682   -0.1598    -0.0583    0.0421   0.0538   -0.0406    0.0038   -0.0118]';
% 
% clearvars explor_X explor_WD_multip
% 
% % FV = [IVs(:,1) week_IV]*[3.0426   12.7745]'; 
% % FV = [IVs(:,1) week_IV innov_X explor_X innov_WD_multip explor_WD_multip]*bs;
% prev_FV_WD = prev_FV_basic + prev_FV_innov + prev_FV_explor;
% 
% save('prev_FV_WD.mat','prev_FV_WD','-v7.3');
% clearvars prev_FV_basic  prev_FV_innov  prev_FV_explor
% 
% display('prev_FV_WD')

% %%
% X = matp(:,[6 7]);
% y = matp(:,5);
% 
% IVs = X;
% choice_dv = [y 1-y];
% 
% week_IV = 100*gampdf(IVs(:,2),0.8375,27.4694);              % w: NOTICE new lines here for fixed gamma parameters
% week_IV(IVs(:,2)<1)=0; 
% 
% clearvars matp
% load('innov_contin_std.mat');
% load('explor_contin_std.mat');
% 
% prev_innov_IV = innov_contin(:,1:4);
% prev_explor_IV = explor_contin(:,1:4);
% clearvars innov_contin explor_contin
% 
% prev_bs = [0.3870    0.1020    -0.0285    0.0031   -0.2426    0.2146    0.1718   -0.1603   -0.0624    0.0459]';
% prev_FV = [IVs(:,1) week_IV prev_innov_IV prev_explor_IV]*prev_bs;
% save('prev_FV.mat','prev_FV','-v7.3');
% clearvars prev_innov_IV prev_explor_IV
% 
% display('prev_FV')

%%
load('matpstd.mat');
% load('innov_contin_std.mat');
% load('explor_contin_std.mat');

X = matp(:,[6 7]);
y = matp(:,5);

IVs = X;
choice_dv = [y 1-y];

clearvars matp

week_IV = 100*gampdf(IVs(:,2),0.8375,27.4694);              % w: NOTICE new lines here for fixed gamma parameters
week_IV(IVs(:,2)<1)=0; 

rb = -10;
ev_points = [];
while rb < 0
    ev_points = [ev_points rb];
    step = abs(rb/10);
    if step < 0.01
        step = 0.01;
    end
    rb = rb + step;
end
ev_points = sort([ev_points -ev_points])
% scatter(1:length(ev_points), ev_points)

seq = [5 6 5 6];
ea_switch = [1 1 0 0];
% prev_joint_bs = [-0.0285    0.0031;   -0.2426    0.2146;    0.1718   -0.1603;   -0.0624    0.0459];

for i = 1:4
    load('innov_contin_std.mat');
    load('explor_contin_std.mat');
    if ea_switch(i) ==1
        ea_cols = seq(i);
        explor_cols = [];
        innov_IV = innov_contin(:,ea_cols);
        explor_IV = [];
    else
        ea_cols = [];
        explor_cols = seq(i);
        innov_IV = [];
        explor_IV = explor_contin(:,explor_cols);
    end   
    
%     ea_cols = [1];
%     % ea_cols = 3;
%     explor_cols = [];
%     innov_IV = innov_contin(:,[1 3]);
%     innov_IV = innov_contin(:,ea_cols)*(prev_joint_bs(i,:)');
%     innov_IV = [];
%     explor_IV = explor_contin(:,explor_cols)*(prev_joint_bs(i,:)');
%     explor_IV = [];
    
    clearvars innov_contin explor_contin
    
%     innov_WD_multip = zeros(size(innov_IV));
%     for i = 1:size(innov_IV,2);
%         innov_WD_multip(:,i) = innov_IV(:,i).*week_IV;
%     end
    innov_WD_multip = [];
    
%     explor_WD_multip = zeros(size(explor_IV));
%     for j = 1:size(explor_IV,2);
%         explor_WD_multip(:,j) = explor_IV(:,j).*week_IV;
%     end
    explor_WD_multip = [];
    
    load('prev_FV_WD.mat');
    
    b_i = ev_points;
    display('LL')
    
    ll_line = zeros(1,length(b_i));
    for d = 1:length(b_i)
        
        const = -5.6271;
        %     bs = [0.3900 0.1015 -0.2641 -0.5259 b_i(d)]';
        %     FV = [IVs(:,1) week_IV innov_IV  innov_WD_multip explor_IV explor_WD_multip]*bs;
        bs = b_i(d);
%         FV = [innov_WD_multip explor_WD_multip].*bs;
        FV = [innov_IV explor_IV].*bs;
        FV = FV + prev_FV_WD;
        
        exp_util = exp(-(const+FV));         % this is now the utility of the external good
        prob=1./(1+exp_util);                % this is still the probability of choosing the product
        pmat = [prob 1-prob];
        pmat = pmat.*choice_dv;
        [r c p] = find(pmat);                                             % I*1
        ll_line(d) = sum(log(p));
    end
    
    
    mat_1d_name = ['ll_line_' num2str(ea_cols) '_' num2str(explor_cols) '_' num2str(not(isempty(innov_WD_multip))) '_' num2str(not(isempty(explor_WD_multip))) '_' num2str(b_i(1)) '_' num2str(b_i(end)) '.mat'];
    save(mat_1d_name,'ll_line') ;
    
    [mv,id] = max(ll_line);
    parameter = b_i(id)
    % [i,j,d] = ind2sub(size(ll_cube),id)
    
    plot(b_i,ll_line)
    
end


% b_i = -10:5:10;
% b_j = -10:5:10;
% display('LL_mat')
% ll_mat = zeros(length(b_i), length(b_j));
% for d = 1:length(b_i)
%     for k = 1:length(b_j)
%         const = -6.1668;
%         bs = [3.0426 12.7745 b_i(d) b_j(k)]';
%         FV = [IVs(:,1) week_IV innov_IV  innov_WD_multip explor_IV explor_WD_multip]*bs;
%         
%         exp_util = exp(-(const+FV));         % this is now the utility of the external good
%         prob=1./(1+exp_util);                % this is still the probability of choosing the product
%         pmat = [prob 1-prob];
%         pmat = pmat.*choice_dv;
%         [r c p] = find(pmat);                                             % I*1
%         ll_mat(d,k) = sum(log(p));
%     end
% end
% 
% mat_1d_name = ['ll_mat_' num2str(ea_cols) '_' num2str(explor_cols) '_' num2str(not(isempty(innov_WD_multip))) '_' num2str(not(isempty(explor_WD_multip))) '_' num2str(b_i(1)) '_' num2str(b_i(end)) '.mat'];
% save(mat_1d_name,'ll_mat') ;
% 
% % [mv,id] = max(ll_line);
% % parameter = b_i(id)
% % [i,j,d] = ind2sub(size(ll_cube),id)
% 
% % [mv,id] = max(ll_mat(:));
% % [i,j,d] = ind2sub(size(ll_cube),id)
% [mv,id] = max(ll_mat(:));
% [i,j] = ind2sub(size(ll_mat),id);
% parameters = [b_i(i) b_i(j)]
% 
% surf(ll_mat)


% %% evaluation points
% X = matp(:,[6 7]);
% y = matp(:,5);
% 
% IVs = X;
% choice_dv = [y 1-y];
% 
% clearvars matp
% 
% week_IV = gampdf(IVs(:,2),0.8343,27.4712);              % w: NOTICE new lines here for fixed gamma parameters
% week_IV(IVs(:,2)<1)=0; 
% 
% rb = -10;
% ev_points = [];
% while rb < 0
%     ev_points = [ev_points rb];
%     step = abs(rb/4);
%     if step < 0.01
%         step = 0.01;
%     end
%     rb = rb + step;
% end
% ev_points = sort([ev_points -ev_points])
% % scatter(1:length(ev_points), ev_points)
% 
% seq = [1 2; 3 4; 1 2; 3 4];
% ea_switch = [1 1 0 0];
% 
% for i = 1:4
%     load('innov_contin_std.mat');
%     load('explor_contin_std.mat');
%     if ea_switch(i) ==1
%         ea_cols = seq(i,:);
%         explor_cols = [];
%     else
%         ea_cols = [];
%         explor_cols = seq(i,:);
%     end
%     innov_IV = innov_contin(:,ea_cols);
%     explor_IV = explor_contin(:,explor_cols);
%     clearvars innov_contin explor_contin
%         
%     % innov_WD_multip = zeros(size(innov_IV));
%     % for i = 1:size(innov_IV,2);
%     %     innov_WD_multip(:,i) = innov_IV(:,i).*week_IV;
%     % end
%     innov_WD_multip = [];
%     
%     % explor_WD_multip = zeros(size(explor_IV));
%     % for j = 1:size(explor_IV,2);
%     %     explor_WD_multip(:,j) = explor_IV(:,j).*week_IV;
%     % end
%     explor_WD_multip = []; 
%     
%     display('LL_mat')
%     
%     ll_mat = zeros(length(ev_points), length(ev_points));
%     for d = 1:length(ev_points)
%         for k = 1:length(ev_points)
%             const = -5.6204;
%             bs = [0.3900 0.1015 ev_points(d) ev_points(k)]';
%             FV = [IVs(:,1) week_IV innov_IV  innov_WD_multip explor_IV explor_WD_multip]*bs;
%             
%             exp_util = exp(-(const+FV));         % this is now the utility of the external good
%             prob=1./(1+exp_util);                % this is still the probability of choosing the product
%             pmat = [prob 1-prob];
%             pmat = pmat.*choice_dv;
%             [r c p] = find(pmat);                                             % I*1
%             ll_mat(d,k) = sum(log(p));
%         end
%     end    
%     mat_1d_name = ['ll_mat_' num2str(ea_cols) '_' num2str(explor_cols) '_' num2str(not(isempty(innov_WD_multip))) '_' num2str(not(isempty(explor_WD_multip))) '_' num2str(ev_points(1)) '_' num2str(length(ev_points)) '_' num2str(ev_points(end)) '.mat'];
%     save(mat_1d_name,'ll_mat') ;
%     
%     [mv,id] = max(ll_mat(:));
%     [i,j] = ind2sub(size(ll_mat),id);
%     parameters = [ev_points(i) ev_points(j)]
%     
%     surf(ll_mat)
% end




% plus = zeros(length(choice_dv),1);
% plus(IVs(:,2)==0)=0.01;               % w: change all 0s to 0.01s, LL increases (for gamma parameter k > 1) ? 
% IVs(:,2) = IVs(:,2)+plus;
% clearvars plus matp

% dummy_X = [member_dummies member_dummies_week_d];
% X = mat1(:,[5 7]);
% y = mat1(:,4);
% beta_0 = zeros(1,size(X,2)+size(dummy_X,2)+1);
% beta_0 = zeros(1,size(X,2)+1);
% beta_0 = zeros(1,size(X,2)+2);

% b_k = 0:0.5:5;
% b_theta = 0:1:10;
% b_wd = 0:0.2:2;     

% b_k = 0.01:0.5:2.01;
% b_theta = 0.01:1:10.01;
% b_wd = 0:1:4;                   % 0.2495

% b_k = 0.01:0.5:4.01;
% b_theta = 0.01:1:7.01;
% b_wd = 2;    
% 
% display('LL')
% 
% ll_cube = zeros(length(b_k),length(b_theta),length(b_wd));
% for i = 1:length(b_k)
%     for j = 1:length(b_theta)
%         for d = 1:length(b_wd)
%             
%             b = [-6.1474 b_k(i) b_theta(j) 3.1975 b_wd(d)];
%             
%             const = b(1);
%             bs = b(4:end)';
%             week_IV = IVs(:,3).*gampdf(IVs(:,2),b_k(i),b_theta(j));
%             FV = [IVs(:,1) week_IV]*bs;
%             
%             exp_util = exp(-(const+FV));         % this is now the utility of the external good
%             prob=1./(1+exp_util);                % this is still the probability of choosing the product
%             pmat = [prob 1-prob];
%             pmat = pmat.*choice_dv;
%             [r c p] = find(pmat);                                             % I*1
%             ll_cube(i,j,d) = sum(log(p));
%         end
%     end
% end
% 
% mat_3d_name = ['ll_cube_' num2str(b_k(1)) '_' num2str(b_k(end)) '_' num2str(b_theta(1)) '_' num2str(b_theta(end)) '_' num2str(b_wd(1)) '_' num2str(b_wd(end)) '.mat'];
% save(mat_3d_name,'ll_cube') ;
% 
% [mv,id] = max(ll_cube(:));
% [i,j,d] = ind2sub(size(ll_cube),id)
% 
% surf(ll_cube(:,:,d))
% % surf(squeeze(ll_cube(i,:,:)))
% % surf(squeeze(ll_cube(:,j,:)))
% 
% parameters = [b_k(i) b_theta(j) b_wd(d)]


% k = find(ll_cube>=-1236400);
% [i,j,d] = ind2sub(size(ll_cube),k);
% ll_cube(i,j,d)=-2000000;
% [mv,id] = max(ll_cube(:));
% [i,j,d] = ind2sub(size(ll_cube),id)
% ll_cube(i,j,d)

% beta_0 = [0 1 2 0]
% beta_0 = [-100 100 -10]
% beta_0 = [-5.6814 2.6098 -0.0040]
% [b, hessian, grad, standard_error, covariance_matrix, t_stat, exit_flag, output] = band_runbi_ll_p(X, y, dummy_X, beta_0);
% [b, hessian, grad, standard_error, covariance_matrix, t_stat, exit_flag, output] = band_runbi_ll_p(X, y, beta_0);
% 
% save('b.mat','b') ;
% save('standard_error.mat','standard_error') ;
% save('t_stat.mat','t_stat') ;
% save('exit_flag.mat','exit_flag') ;
% 
% display(b)
% % display(standard_error)
% display(t_stat)
% display(grad)
% display(output)
% 
% m = grad'*(-inv(hessian))*grad;          % convergence criterion of Train

% diary off 

%%
% load('row_mid.mat');
% load('row_num.mat');
% 
% I = 158
% 
% b_store_p = cell(I,1);
% se_store_p = cell(I,1);
% cov_store_p = cell(I,1);
% tstat_store_p = cell(I,1);
% % hess = cell(I,1);
% i_id = zeros(I,1);
% est_exitflag_p = zeros(I,1);
% 
% ind = 0;
% for i = 1:I
%     i_id(i) = row_mid(i);
%     nrows = row_num(i);
%     imat = matp(ind+1:ind+nrows,:);         % ind tbc
%     
%     X = imat(:,[6 7]);
%     y = imat(:,5);
%     beta_0 = zeros(1,size(X,2)+1);
%     [b, hessian, standard_error, covariance_matrix, t_stat, exit_flag] = band_runbi_ll_p(X, y, beta_0);
% 
%     b_store_p{i} = b;
%     se_store_p{i} = standard_error';
%     cov_store_p{i} = covariance_matrix;
%     tstat_store_p{i} = t_stat;
% %    hess{i} = hessian;
%     est_exitflag_p(i) = exit_flag;
%         
%     ind = ind+nrows; 
% %     if i >= 10
% %         break
% %     end    
% end
%    
% save('b_store_p.mat','b_store_p') ;
% save('se_store_p.mat','se_store_p') ;
% save('tstat_store_p.mat','tstat_store_p') ;
% save('est_exitflag_p.mat','est_exitflag_p') ;

% b_k = 0:0.5:5;
% b_theta = 0:1:10;
% b_wd = 0:0.2:4;  
% 
% load('ll_cube_0.01_2.01_0.01_10.01_0_4.mat')
% ll_cube3 = ll_cube;
% 
% load('ll_cube_0_5_0_10_0_2.mat')
% ll_cube2 = ll_cube;
% 
% N = NaN(size(ll_cube2,1),size(ll_cube2,2),size(ll_cube2,3));
% ll_combine = cat(3,ll_cube2,N);
% 
% for i = 1:5
%     ll_combine(1:5,:,(i-1)*5+1) = ll_cube3(:,:,i);
% end
% 
% [mv,id] = max(ll_combine(:));
% [i,j,d] = ind2sub(size(ll_combine),id)
% 
% surf(ll_combine(:,:,d))
% parameters = [b_k(i) b_theta(j) b_wd(d)]
% 
% % % no need for this, can just use A(A>0.5)=1
% % k = find(ll_combine<-1300000);
% % [i,j,d] = ind2sub(size(ll_combine),k);
% % for ind = 1:length(k)
% %     ll_combine(i(ind),j(ind),d(ind))=NaN;
% % end
% ll_combine(ll_combine<-1250000)=NaN;
% [mv,id] = max(ll_combine(:));
% [i,j,d] = ind2sub(size(ll_combine),id)
% ll_combine(i,j,d)
% surf(ll_combine(:,:,d))
%     
% % surf(ll_combine(3:11,2:11,11))
