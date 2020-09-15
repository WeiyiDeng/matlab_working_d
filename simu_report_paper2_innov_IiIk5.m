clc
clear

rng(181102)

% date_=clock;
% resultsfilename=['Results/r_Results_' num2str(date_(1)) '_' num2str(date_(2)) '_' num2str(date_(3)) '-' num2str(date_(4))  '_' num2str(date_(5)) '.txt'];
% diary(resultsfilename);
% 
% mfilename
% p = mfilename('fullpath')

load('matp_friend_reverse.mat');
predict_trend = csvread('predict_trend_log4061_lenient.csv',1,0);

index = find(ismember(matp(:,3),predict_trend(:,1)));
matp = matp(index,:);

% prep_matp = matp(:,[1 3 4]);
% 
% dummy_innov_mat = csvread('dummy_SI_innov_mat.csv');
% temp = dummy_innov_mat(end,:);
% dummy_innov_mat = [temp; dummy_innov_mat];
% dummy_innov_mat(end,:) = [];
% dum_mat_innov_prep = dummy_innov_mat(:,1:3);
% 
% [~,indx]=ismember(dum_mat_innov_prep,prep_matp,'rows');
% 
% save('indx_innov.mat','indx') ;
% save('dummy_innov_mat.mat','dummy_innov_mat') ;

%
load('indx.mat');
load('dummy_mat.mat');

matp = matp(indx,:);

X = dummy_mat(:, 5);
y = dummy_mat(:,4);
dummy_agg_SI = dummy_mat(:,6);

matp_m_id = matp(:,1);

load('mat_sparse_8088.mat');
load('m_innov.mat');
load('f_innov.mat');
load('cosine_similarity_scores_friends8088_listen1.mat');

hist(f_innov)
hist(m_innov)
mean(m_innov)
mean(f_innov)
quantile(m_innov,0.1)
quantile(m_innov,0.9)
quantile(m_innov,0.3)
quantile(f_innov,0.9)
quantile(f_innov,0.3)

% rm_lowInnovf = f_innov > quantile(f_innov,0.3);
% rm_highInnovf = f_innov < quantile(f_innov,0.6);
% rm_mediumInnovf = f_innov > quantile(f_innov,0.6) | f_innov < quantile(f_innov,0.3);

friendlist8088 = csvread('new_friendlist_8088.csv',1,0);
friendlist8088 = [friendlist8088(end,:); friendlist8088];
friendlist8088(end,:) = [];

m_uni = unique(friendlist8088(:,1));
innov_m_uni = zeros(size(m_uni));
for i = 1:length(m_uni)
    inds = find(friendlist8088(:,1)==m_uni(i));
    innov_m_uni(i) = m_innov(inds(1));
end
% 
% lowInnovm = innov_m_uni > quantile(innov_m_uni,0.3);
% highInnovm = innov_m_uni < quantile(innov_m_uni,0.6);
% mediumInnovm = innov_m_uni > quantile(innov_m_uni,0.6) | innov_m_uni < quantile(innov_m_uni,0.3);

temp = unique(friendlist8088(:));
innov_uni = zeros(size(temp));
for i = 1:length(temp)
    inds = find(friendlist8088(:)==temp(i));
    merge_innov = [m_innov; f_innov];
    innov_uni(i) = merge_innov(inds(1));
end

rm_lowInnovf = f_innov > quantile(innov_uni,0.33);
rm_highInnovf = f_innov < quantile(innov_uni,0.66);
rm_mediumInnovf = f_innov > quantile(innov_uni,0.66) | f_innov < quantile(innov_uni,0.33);

lowInnovm_ind = find(innov_m_uni < quantile(innov_uni,0.33));
highInnovm_ind = find(innov_m_uni > quantile(innov_uni,0.66));
mediumInnovm_ind = find(innov_m_uni < quantile(innov_uni,0.66) & innov_m_uni > quantile(innov_uni,0.33));
lowInnovm = m_uni(lowInnovm_ind);
highInnovm = m_uni(highInnovm_ind);
mediumInnovm = m_uni(mediumInnovm_ind);
lowInnovm_ind_long = ismember(matp_m_id, lowInnovm);
highInnovm_ind_long = ismember(matp_m_id, highInnovm);
mediumInnovm_ind_long = ismember(matp_m_id, mediumInnovm);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rm_choice_Innovf = rm_lowInnovf;
% rm_choice_Innovf = rm_mediumInnovf;
% rm_choice_Innovf = rm_highInnovf;
rm_choice_Innovf = ones(size(f_innov));

load('Im_17617085.mat');
load('Ik_17617085.mat');

% X = dummy_innov_mat(:, 5);
% y = dummy_innov_mat(:,4);
% dummy_agg_SI = dummy_innov_mat(:, 6);
% dummy_agg_SI_innov = dummy_innov_mat(:, 8);

% load('PI_member.mat');

% temp_memberPI = sparse(1,PI_member(:,1),PI_member(:,2),1,8320);

% memberPI_vec = zeros(size(matp(:,1),1),1);
% for i = 1:length(memberPI_vec)
%     memberPI_vec(i) = temp_memberPI(matp(i,1));
% end

%%
band_pop = csvread('POP_counts_year_final.csv',1,0);
band_pop(:,2) = band_pop(:,2)-104;
band_pop_mat = sparse(band_pop(:,1),band_pop(:,2),band_pop(:,4),max(band_pop(:,1)),max(band_pop(:,2)));
pop = zeros(size(matp,1),1);
for r = 1:length(pop)
    pop(r) = band_pop_mat(matp(r,3),matp(r,4));
end

predict_trend_mat = sparse(predict_trend(:,1),predict_trend(:,2),predict_trend(:,3),...
    max(predict_trend(:,1)),max(predict_trend(:,2)));

trend_hat = zeros(size(matp,1),1);
for r = 1:length(trend_hat)
    trend_hat(r) = predict_trend_mat(matp(r,3),matp(r,4));
end
matp(:,6) = trend_hat;

corr(trend_hat,pop)

band_topics = csvread('band_count_topics15.csv',1,0);
band_topics_mat = sparse(band_topics(:,1),1,band_topics(:,2),max(band_topics(:,1)),1);
topics_count = zeros(size(matp,1),1);
for r = 1:length(topics_count)
    topics_count(r) = band_topics_mat(matp(r,3));
end

band_birth = csvread('introdate3.csv',1,0);
band_birth_mat = sparse(band_birth(:,1),1,band_birth(:,2),max(band_birth(:,1)),1);
band_age = zeros(size(matp,1),1);
for r = 1:length(band_age)
    band_age(r) = matp(r,4)-band_birth_mat(matp(r,3));
end
sum(band_age<0)

band_count_tracks = csvread('BAND_count_TRACKS.csv',1,0);
band_tracks_mat = sparse(band_count_tracks(:,1),1,band_count_tracks(:,2),max(band_count_tracks(:,1)),1);
tracks_count = zeros(size(matp,1),1);
for r = 1:length(tracks_count)
    tracks_count(r) = matp(r,4)-band_tracks_mat(matp(r,3));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('N_prep_for_matp_jobs_organize.mat');
load('combi_similarities.mat')

N_prep_for_matp_jobs_organize = N_prep_for_matp_jobs_organize(index,:);
N_prep_for_matp_jobs_organize = N_prep_for_matp_jobs_organize(indx,:);

[row col val] = find(N_prep_for_matp_jobs_organize==1);
store_r = zeros(length(row),1);
for i = 1:length(row)-1
    if row(i+1) == row(i)+2
        store_r(i) = 1;
    end
end
row0 = row(store_r==1)+1;
col0 = col(store_r==1);
val0 = zeros(size(col0));

row_store = [row; row0];
col_store = [col; col0];
val_store = [val; val0];

% val_pdf = normpdf(val_store);
% tic
% N_prep_for_matp_jobs_organize=sparse(row_store,col_store,val_pdf);
% toc

None0s_X_N = [row_store col_store val_store];
X_N = N_prep_for_matp_jobs_organize;
S = combi_similarities;

clearvars N_prep_for_matp_jobs_organize combi_similarities matp

%%
load('b_agg_dummy_pop3_full.mat')
% beta_0 = b
% beta_0 = [-7.7573    0.0030   -0.0004    0.0152   -0.0282    0.0230    0.0139    0.0489    2.1851    0.7553]
% beta_0 = [-7.7250    0    0    0.01    0.01    0.01    0.01    0    2.1842   0.7598   0]
% beta_0 = [-7.7250    0    0    0.01    0.01    0.01    0.01    0    2.1842   0.7598]
% beta_0 = [-7.7250    0.01    0.01    0.01    0    0    2.1842   0.7598]
% clearvars b

const = b(1);

% FV = IVs*bs;
b_basic = b(2:8)';

% week_IV = dummy_agg_SI;
% week_IV_innov = dummy_agg_SI_innov;

Aijkt = mat_sparse_8088;
Sik = cosine_similarity_scores;

val_pdf = 100*normpdf(None0s_X_N(:,3),0,b(9));
X_N = sparse(None0s_X_N(:,1),None0s_X_N(:,2),val_pdf,17617085,6222);
IV_N_S = X_N*S.^exp(b(10));

pop = pop./1000;

f_innov = f_innov*10;
m_innov = m_innov*10;
Aijkt = Aijkt*10;

% Sik_power = Sik.^exp(b(11));
% Sik_power = Sik;
Sik_power = ones(8088,1);
% IV_N_S = IV_N_S*10;

% Im_17617085 = Im_17617085./10;
% Ik_17617085 = Ik_17617085./10;

% FV = [IVs(:,1)  trend_hat  week_IV  band_age  topics_count...
%     band_age.*week_IV  topics_count.*week_IV   band_age.*topics_count...
%     band_age.*topics_count.*week_IV    IV_N_S]*b_basic;          % with both trend and BP as controls
% FV = [trend_hat  week_IV  week_IV.^2   pop   week_IV.*pop      week_IV.^2.*pop     week_IV_innov    IV_N_S]*b_basic;          % with both trend and BP as controls
FV = [Im_17617085    Ik_17617085     Aijkt*Sik_power    Aijkt*(m_innov.*Sik_power)    Aijkt*(f_innov.*Sik_power.*rm_choice_Innovf)    Aijkt*(f_innov.*m_innov.*Sik_power.*rm_choice_Innovf)     IV_N_S]*b_basic;

% week_IV.^2*pop

% exp_util = exp(const+FV);          % utility of choosing the product
% prob=exp_util./(1+exp_util);
exp_util = exp(-(const+FV));         % this is now the utility of the external good
prob=1./(1+exp_util);                % this is still the probability of choosing the product
% pmat = [prob 1-prob]; 

iter = 300;
num_sample = zeros(1,iter);
num_Lowinnovm = zeros(1,iter);
num_Mediuminnovm = zeros(1,iter);
num_Highinnovm = zeros(1,iter);
for i = 1:iter
    temp = rand(size(y));
    num_sample(i) = sum(temp<prob);
    num_Lowinnovm(i) = sum((temp<prob).*lowInnovm_ind_long);
    num_Mediuminnovm(i) = sum((temp<prob).*mediumInnovm_ind_long);
    num_Highinnovm(i) = sum((temp<prob).*highInnovm_ind_long);
end
mean(num_sample)
mean(num_Lowinnovm)
mean(num_Mediuminnovm)
mean(num_Highinnovm)

% confusion matrix
pred = temp<prob;
sum(pred)
sum(pred==y & y==1)

confusion_matrix = [sum(pred==y & y==1) sum(pred~=y & y==0); sum(pred~=y & y==1) sum(pred==y & y==0)]
accuracy = (sum(pred==y & y==1)+sum(pred==y & y==0))/sum(sum(confusion_matrix))

simu_sample_all_lowinnovm_mediuminovm_highinnovm = [mean(num_sample) mean(num_Lowinnovm) mean(num_Mediuminnovm) mean(num_Highinnovm)];
% save('simu_sample_rm_lowInnovf.mat','simu_sample_all_lowinnovm_mediuminovm_highinnovm');
% save('simu_sample_rm_mediumInnovf.mat','simu_sample_all_lowinnovm_mediuminovm_highinnovm');
% save('simu_sample_rm_highInnovf.mat','simu_sample_all_lowinnovm_mediuminovm_highinnovm');
% save('simu_sample_rm_geenInnovf.mat','simu_sample_all_lowinnovm_mediuminovm_highinnovm');

% 8.1211e+03    rm_H
% 1.0334e+04    rm_M
% 1.2782e+04    rm_L
% 1.3474e+04    non_rm

load('simu_sample_rm_lowInnovf.mat')
load('simu_sample_rm_mediumInnovf.mat')
load('simu_sample_rm_highInnovf.mat')
load('simu_sample_rm_geenInnovf.mat')

% [b, hessian, grad, standard_error, covariance_matrix, t_stat, exit_flag, output] = simu_runbi_ll_paper2_innov_IiIk5(X, trend_hat, pop, dummy_agg_SI, mat_sparse_8088, m_innov, f_innov, Im_17617085, Ik_17617085, cosine_similarity_scores, None0s_X_N, S, y, beta_0, rm_choice_Innovf);
% 
% save('b_agg_dummy_pop3_full.mat','b') ;
% save('standard_error_agg_dummy_pop3_full.mat','standard_error') ;
% save('t_stat_agg_dummy_pop3_full.mat','t_stat') ;
% save('exit_flag_agg_dummy_pop3_full.mat','exit_flag') ;
% save('hessian_pop3_full.mat','hessian') ;
% 
% display(b)
% % display(standard_error)
% display(t_stat)
% display(grad)
% display(output)

% diary off





