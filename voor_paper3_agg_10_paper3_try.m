clc
clear

load('matp_friend_reverse.mat');
predict_trend = csvread('predict_trend_log4061_lenient.csv',1,0);

index = find(ismember(matp(:,3),predict_trend(:,1)));

matp = matp(index,:);

% load('dummies40.mat');
% % size(dummies40)
% dummies40 = dummies40(index,:);
% dummies_prep = dummies40.*repmat(matp(:,8),1,4);
% clearvars dummies40
% dummy_prep = dummies_prep(:,1) | dummies_prep(:,2);         % 1-20
% 
% % mat_export = [matp(:,1:6) dummy_prep];
% % csvwrite('mat_export.csv',mat_export)
% 
% prep_matp = matp(:,[1 3 4]);
% 
% % remove duplicte commas in EMeditor
% dummy_mat = csvread('dummy_SI_mat.csv');
% dum_mat_prep = dummy_mat(:,1:3);
% 
% [~,indx]=ismember(dum_mat_prep,prep_matp,'rows');
% 
% save('indx.mat','indx') ;
% save('dummy_mat.mat','dummy_mat') ;
load('indx.mat');
load('dummy_mat_fix.mat');

matp = matp(indx,:);

% X = dummy_mat(:, 5);
DV_paper1 = dummy_mat(:,4);
dummy_agg_SI = dummy_mat(:, 6);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

matp_mjt = matp(:,[1 3 4]);

clearvars N_prep_for_matp_jobs_organize combi_similarities matp

%%
load('b_agg_dummy_pop3_full_fix.mat')
% beta_0 = b
% b = [-8.0906    0.3280    0.2663   -0.0070    0.8637    -0.1342    0.0048   0.0391    2.1956   0.7717]
% clearvars b

b_basic = b(2:8)';

week_IV = dummy_agg_SI;

val_pdf = 100*normpdf(None0s_X_N(:,3),0,b(9));
X_N = sparse(None0s_X_N(:,1),None0s_X_N(:,2),val_pdf,17617085,6222);
IV_N_S = X_N*S.^exp(b(10));

pop = pop./1000;

% FV = [IVs(:,1)  trend_hat  week_IV  band_age  topics_count...
%     band_age.*week_IV  topics_count.*week_IV   band_age.*topics_count...
%     band_age.*topics_count.*week_IV    IV_N_S]*b_basic;          % with both trend and BP as controls
FV = [trend_hat  week_IV  week_IV.^2   pop   week_IV.*pop      week_IV.^2.*pop     IV_N_S]*b_basic;

prep_SI_pop_paper3 = [week_IV.*b(3)  week_IV.^2.*b(4)   pop.*b(5)   week_IV.*pop.*b(6)      week_IV.^2.*pop.*b(7)];

prep_paper3 = [matp_mjt prep_SI_pop_paper3];
prep_paper3_fix = [matp_mjt prep_SI_pop_paper3 DV_paper1]; 
% save('prep_paper3.mat','prep_paper3') 
% csvwrite('prep_paper3.csv',prep_paper3)
% csvwrite('prep_paper3_fix.csv',prep_paper3_fix)

%%

m_t_agg = csvread('agg10_paper3_mat.csv');
m_t_usage = csvread('dt7_member_ID_sum_listen.csv',1,0);

m_t_usage(:,2) = m_t_usage(:,2)-104;
ind0 = m_t_usage(:,2)<=0;
m_t_usage(ind0,:) = [];
usage_sparse = sparse(m_t_usage(:,1),m_t_usage(:,2),m_t_usage(:,3),8320,423);
ind423 = m_t_agg(:,2)==423;
m_t_agg(ind423,:) = [];
paper3_DV = zeros(size(m_t_agg,1),1);
paper3_DV_lag = zeros(size(m_t_agg,1),1);
for i = 1:size(m_t_agg,1)
    paper3_DV(i) = usage_sparse(m_t_agg(i,1),m_t_agg(i,2)+1);
    paper3_DV_lag(i) = usage_sparse(m_t_agg(i,1),m_t_agg(i,2));
end

% for i = 1:size(m_t_agg,1)
%     ind = find(m_t_usage(:,1)==m_t_agg(i,1) & m_t_usage(:,2)==m_t_agg(i,2));
%     if length(ind)>1
%         disp(i)
%     elseif isempty(ind)
%         paper3_DV(i) = 0;
%     else
%         paper3_DV(i) = m_t_usage(ind,3);
%     end
% end

%%
load('indx.mat');
load('dummy_mat_fix.mat');

% X = dummy_mat(:, 5);
DV_paper1 = dummy_mat(:,4);
dummy_agg_SI = dummy_mat(:, 6);

collect_dummy_mat = dummy_mat(find(dummy_mat(:,4)),:);
sparse_dummy = sparse(collect_dummy_mat(:,1),collect_dummy_mat(:,3),1,8320,423);
sum_paper1_DV_dummy = zeros(size(m_t_agg,1),1);
for i = 1:size(m_t_agg,1)
    sum_paper1_DV_dummy(i) = sparse_dummy(m_t_agg(i,1),m_t_agg(i,2));
end

%%
sum_paper1_DV_dummy_ASI_APOP = sum_paper1_DV_dummy - m_t_agg(:,3)-m_t_agg(:,4);
% X = [m_t_agg(:,3:6) paper3_DV_lag sum_paper1_DV_dummy_ASI_APOP];
% X = [m_t_agg(:,3:4) paper3_DV_lag];
% X = [m_t_agg(:,5:6) paper3_DV_lag sum_paper1_DV_dummy_ASI_APOP];
% X = [m_t_agg(:,3:4) paper3_DV_lag];
% X = [m_t_agg(:,3) m_t_agg(:,4) sum_paper1_DV_dummy m_t_agg(:,3).*sum_paper1_DV_dummy m_t_agg(:,4).*sum_paper1_DV_dummy];
X = [m_t_agg(:,3) m_t_agg(:,4) paper3_DV_lag m_t_agg(:,3).*paper3_DV_lag m_t_agg(:,4).*paper3_DV_lag];
% X = [m_t_agg(:,3) m_t_agg(:,4) paper3_DV_lag sum_paper1_DV_dummy m_t_agg(:,3).*sum_paper1_DV_dummy m_t_agg(:,4).*sum_paper1_DV_dummy];

X_std = (X-mean(X))./std(X);
y = paper3_DV;

corr(X)
mdl = fitlm(X_std,y)

csvwrite('X_std_paper3_Nov.csv',X_std)
csvwrite('y_paper3_Nov.csv',y)

%%
load('m_innov.mat');



