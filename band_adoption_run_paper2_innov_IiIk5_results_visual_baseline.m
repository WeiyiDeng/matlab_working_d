clc
clear

% date_=clock;
% resultsfilename=['Results/r_Results_' num2str(date_(1)) '_' num2str(date_(2)) '_' num2str(date_(3)) '-' num2str(date_(4))  '_' num2str(date_(5)) '.txt'];
% diary(resultsfilename);

mfilename
p = mfilename('fullpath')

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

choice_dv = [y 1-y];

load('mat_sparse_8088.mat');
load('m_innov.mat');
load('f_innov.mat');
load('cosine_similarity_scores_friends8088_listen1.mat');

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
% load('b_agg_dummy_pop.mat')
% beta_0 = b
% beta_0 = [-7.7573    0.0030   -0.0004    0.0152   -0.0282    0.0230    0.0139    0.0489    2.1851    0.7553]
% beta_0 = [-7.7250    0    0    0.01    0.01    0.01    0.01    0    2.1842   0.7598   0]
% beta_0 = [-7.7250    0    0    0.01    0.01    0.01    0.01    0    2.1842   0.7598]
% beta_0 = [-7.7250    0.01    0.01    0.01    0    0    2.1842   0.7598]
b = [-7.7573    0.0030   -0.0004    0.0152   -0.0282    0.0230    0.0139    0.0489    2.1851    0.7553];

%%
f_innov = f_innov*10;
m_innov = m_innov*10;
% S_percent = quantile(S,[.33 .50 .66]);
% f_innov_percent = quantile(f_innov,[.33 .50 .66])
% m_innov_percent = quantile(m_innov,[.33 .50 .66])
% f_innov_percent = quantile(f_innov,0.1:0.1:1)
% m_innov_percent = quantile(m_innov,0.1:0.1:1)
% f_innov_percent = quantile(f_innov,[.33 .66 .99])
% m_innov_percent = quantile(m_innov,[.33 .66 .99])
f_innov_percent = quantile(f_innov,[.16 .50 .84])
m_innov_percent = quantile(m_innov,[.16 .50 .84])
% f_innov_percent = quantile(f_innov,0.16:0.02:0.84)
% m_innov_percent = quantile(m_innov,0.16:0.02:0.84)

Aijkt = mat_sparse_8088;
% Sik_power = ones(8088,1);
Sik = cosine_similarity_scores;

% Im_17617085 = Im_17617085./10;
% Ik_17617085 = Ik_17617085./10;

% Im_percent = quantile(Im_17617085,[.33 .50 .66])
% Ik_percent = quantile(Ik_17617085,[.33 .50 .66])
Im_percent = f_innov_percent./100;
Ik_percent = m_innov_percent./100;

val_pdf = 100*normpdf(None0s_X_N(:,3),0,b(9));
X_N = sparse(None0s_X_N(:,1),None0s_X_N(:,2),val_pdf,17617085,6222);
IV_N_S = X_N*S.^exp(b(10));
% N_percent = quantile(IV_N_S,[0.33 0.97 0.99])
% N_percent = 0;
% N_percent = mean(IV_N_S)

% Ik_percent = Ik_percent*10;
% Im_percent = Im_percent*10;
Aijkt = Aijkt*10;

% Sik_power = Sik.^exp(b(11)/1000);
Sik_power = Sik;
% SI_percent = quantile(Aijkt*Sik_power,[0.33 0.8 0.99])
% SI_percent = 0;
% SI_percent = mean(Aijkt*Sik_power);

const = b(1);
b_basic = b(2:8)';

temp = Aijkt*Sik_power;
v1 = IV_N_S(find(IV_N_S));
v2 = temp(find(temp));
% N_percent = quantile(v1,[0.33 0.66 0.99]);
% SI_percent = quantile(v2,[0.33 0.66 0.99]);
N_percent = quantile(v1,0.5);                     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SI_percent = quantile(v2,0.5);

ind = 1
bar = length(f_innov_percent);
prob = zeros(bar,bar);
for m = 1:bar
    ind_m = m;
    for k = 1:bar
        ind_k = k;
        % FV = [Im_percent(ind)    Ik_percent(ind)     SI_percent(ind)    SI_percent(ind)*m_innov_percent(ind)    SI_percent(ind)*f_innov_percent(ind)    SI_percent(ind)*m_innov_percent(ind)*f_innov_percent(ind)     N_percent(ind)]*b_basic;
        FV = [Im_percent(ind_m)    Ik_percent(ind_k)     SI_percent(ind)    SI_percent(ind)*m_innov_percent(ind_m)    SI_percent(ind)*f_innov_percent(ind_k)    SI_percent(ind)*m_innov_percent(ind_m)*f_innov_percent(ind_k)     N_percent(ind)]*b_basic;
%         FV = [Im_percent(ind_m)    Ik_percent(ind_k)     0    0    0    0     N_percent(ind)]*b_basic;
        exp_util = exp(-(const+FV));         % this is now the utility of the external good
        prob(m,k)=1./(1+exp_util)                % this is still the probability of choosing the product
    end
end
surf(prob)
xlabel('innov F (sender)')
ylabel('innov M (receiver)')
zlabel('P')
%%
% const = b(1);

% FV = IVs*bs;
% b_basic = b(2:8)';

% week_IV = dummy_agg_SI;
% week_IV_innov = dummy_agg_SI_innov;

% Aijkt = mat_sparse_8088;
% Sik = cosine_similarity_scores;

% val_pdf = 100*normpdf(None0s_X_N(:,3),0,b(9));
% X_N = sparse(None0s_X_N(:,1),None0s_X_N(:,2),val_pdf,17617085,6222);
% IV_N_S = X_N*S.^exp(b(10));

% pop = pop./1000;

% f_innov = f_innov*10;
% m_innov = m_innov*10;
% Aijkt = Aijkt*10;

% Sik_power = Sik.^exp(b(11));
% Sik_power = Sik;
% Sik_power = ones(8088,1);
% IV_N_S = IV_N_S*10;

% Im_17617085 = Im_17617085./10;
% Ik_17617085 = Ik_17617085./10;

% FV = [IVs(:,1)  trend_hat  week_IV  band_age  topics_count...
%     band_age.*week_IV  topics_count.*week_IV   band_age.*topics_count...
%     band_age.*topics_count.*week_IV    IV_N_S]*b_basic;          % with both trend and BP as controls
% FV = [trend_hat  week_IV  week_IV.^2   pop   week_IV.*pop      week_IV.^2.*pop     week_IV_innov    IV_N_S]*b_basic;          % with both trend and BP as controls
% FV = [Im_17617085    Ik_17617085     Aijkt*Sik_power    Aijkt*(m_innov.*Sik_power)    Aijkt*(f_innov.*Sik_power)    Aijkt*(f_innov.*m_innov.*Sik_power)     IV_N_S]*b_basic;

% week_IV.^2*pop

% exp_util = exp(const+FV);          % utility of choosing the product
% prob=exp_util./(1+exp_util);
% exp_util = exp(-(const+FV));         % this is now the utility of the external good
% prob=1./(1+exp_util);                % this is still the probability of choosing the product
% pmat = [prob 1-prob]; 
% pmat = pmat.*choice_dv;
% [r c p] = find(pmat);                                             % I*1
% LL = -sum(log(p));                                                % 1*1

%%
% clearvars b

% [b, hessian, grad, standard_error, covariance_matrix, t_stat, exit_flag, output] = band_runbi_ll_paper2_innov_IiIk5(X, trend_hat, pop, dummy_agg_SI, mat_sparse_8088, m_innov, f_innov, Im_17617085, Ik_17617085, cosine_similarity_scores, None0s_X_N, S, y, beta_0);
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

x1 = prob(1,:);
x2 = prob(2,:);
x3 = prob(3,:);
y = [0.16	0.5	   0.84]
plot(y,x1,'-.o','color','k')
hold on
plot(y,x2,'--o','color','k')
hold on
plot(y,x3,':o','color','k')
hold off
xticks([0.16 0.5 0.84])
legend('Influencee Innovativeness percentile 16%','Influencee Innovativeness percentile 50%','Influencee Innovativeness percentile 84%','Location','northwest')
ylabel('band sampling probability') 
xlabel('Influencer Innovativeness percentile (16%, 50%, 84%)')




