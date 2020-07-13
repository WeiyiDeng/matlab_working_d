clc
clear

rng(181102)

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

temp = matp(17574610:end,:);
matp(17574610:end,:) = [];
matp = [temp; matp];

% X = dummy_mat(:, 5);
DV_paper1 = dummy_mat(:,4);
dummy_agg_SI = dummy_mat(:, 6);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
band_pop = csvread('POP_counts_year_final.csv',1,0);
band_pop(:,2) = band_pop(:,2)-104;

band_pop_adjust_1 = band_pop;
temp_ind = find(band_pop_adjust_1(:,2)==423);
band_pop_adjust_1(temp_ind,:) = [];
band_pop_adjust_1(:,2) = band_pop_adjust_1(:,2)+1;

band_pop_mat = sparse(band_pop(:,1),band_pop(:,2),band_pop(:,4),max(band_pop(:,1)),max(band_pop(:,2)));
pop = zeros(size(matp,1),1);
for r = 1:length(pop)
    pop(r) = band_pop_mat(matp(r,3),matp(r,4));
end

band_pop_mat_adjust_1 = sparse(band_pop_adjust_1(:,1),band_pop_adjust_1(:,2),band_pop_adjust_1(:,4),max(band_pop_adjust_1(:,1)),max(band_pop_adjust_1(:,2)));
pop_adjust_1 = zeros(size(matp,1),1);
for r = 1:length(pop_adjust_1)
    pop_adjust_1(r) = band_pop_mat_adjust_1(matp(r,3),matp(r,4));
end

pop_adjust_1 = pop_adjust_1./1000;

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
const = b(1);

week_IV = dummy_agg_SI;

val_pdf = 100*normpdf(None0s_X_N(:,3),0,b(9));
X_N = sparse(None0s_X_N(:,1),None0s_X_N(:,2),val_pdf,17617085,6222);
IV_N_S = X_N*S.^exp(b(10));

pop = pop./1000;

% FV = [IVs(:,1)  trend_hat  week_IV  band_age  topics_count...
%     band_age.*week_IV  topics_count.*week_IV   band_age.*topics_count...
%     band_age.*topics_count.*week_IV    IV_N_S]*b_basic;          % with both trend and BP as controls
FV = [trend_hat  week_IV  week_IV.^2   pop   week_IV.*pop      week_IV.^2.*pop     IV_N_S]*b_basic;
exp_util = exp(-(const+FV));         % this is now the utility of the external good
prob=1./(1+exp_util);    

%% WHERE THE SIMULATION STARTS
simu_no_SI_FV = [trend_hat  pop      IV_N_S]*b_basic([1 4 7]);
exp_util_noSI = exp(-(const+simu_no_SI_FV));         % this is now the utility of the external good
prob_noSI=1./(1+exp_util_noSI);    

runs = 300
simu_paper3 = cell(runs,1);
for s = 1:runs
    temp = rand(size(prob));
    simu_sample = temp<prob;
    simu_sample_noSI = temp<prob_noSI;
    sum(simu_sample)
    sum(simu_sample_noSI)

    APOP_prep = pop_adjust_1.*simu_sample_noSI;

    max(matp_mjt)
    matp_mjt_rm_ind = find(matp_mjt(:,3)==423);
    matp_mjt_rm423 = matp_mjt;
    matp_mjt_rm423(matp_mjt_rm_ind,:)=[];
    APOP_prep_rm423 = APOP_prep;

    length(unique(matp_mjt(:,1)))

    member_id = matp_mjt_rm423(1,1);
    m_end = [];
    m_start = 1;
    m_list = [];
    for i = 2:size(matp_mjt_rm423,1)
        prev_m = member_id;
        member_id = matp_mjt_rm423(i,1);
        if prev_m ~=member_id
            m_end = [m_end; i-1];
            m_start = [m_start; i];
            m_list = [m_list; prev_m];
        else
        end
    end
    m_list = [m_list; member_id];    
    m_end = [m_end; i]; 

    ind_row = 1:length(m_list);
    member_nummer = sparse(m_list,1,ind_row,8320,1);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    member_APOP_cells = cell(length(m_list),1);
    member_agg_bands_APOP_cells = cell(length(m_list),1);
    % for i = 1:length(unique(matp_mjt(:,1)))
    for i = 1:length(m_list)
        member_APOP_cells{i} = ...
            sparse(matp_mjt_rm423(m_start(i):m_end(i)-1,3)+1,matp_mjt_rm423(m_start(i):m_end(i)-1,2),APOP_prep_rm423(m_start(i):m_end(i)-1),423,6046);
        member_agg_bands_APOP_cells{i} = sum(member_APOP_cells{i},2);
    end

    % % save('simu_sample_noSI.mat','simu_sample_noSI');
    % 
    % prep_SI_pop_paper3 = [week_IV.*b(3)  week_IV.^2.*b(4)   pop.*b(5)   week_IV.*pop.*b(6)      week_IV.^2.*pop.*b(7)];
    % 
    % prep_paper3 = [matp_mjt prep_SI_pop_paper3];
    % prep_paper3_fix = [matp_mjt prep_SI_pop_paper3 DV_paper1]; 
    % % save('prep_paper3.mat','prep_paper3') 
    % % csvwrite('prep_paper3.csv',prep_paper3)
    % % csvwrite('prep_paper3_fix.csv',prep_paper3_fix)

    %%
    % m_t_agg = csvread('agg10_paper3_mat.csv');
    m_t_agg = csvread('agg10_paper3_mat.csv');
    m_t_usage = csvread('dt7_member_ID_sum_listen.csv',1,0);
    load('paper3_stage2_modelresults.mat')

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

    %% SIMU 300 runs
    APOP_ready_vec = zeros(size(m_t_agg,1),1);
    for r = 1:size(m_t_agg,1)
        APOP_ready_vec(r) = member_agg_bands_APOP_cells{member_nummer(m_t_agg(r,1))}(m_t_agg(r,2));
    end

    member_id_rec = m_t_agg(1,1);
    m_end_rec = [];
    m_start_rec = 1;
    m_list_rec = [];
    for i = 2:size(m_t_agg,1)
        prev_m = member_id_rec;
        member_id_rec = m_t_agg(i,1);
        if prev_m ~=member_id_rec
            m_end_rec = [m_end_rec; i-1];
            m_start_rec = [m_start_rec; i];
            m_list_rec = [m_list_rec; prev_m];
        else
        end
    end
    m_list_rec = [m_list_rec; member_id_rec];    
    m_end_rec = [m_end_rec; i];

    min(m_end_rec-m_start_rec)

    SimuX = [ones(size(paper3_DV_lag)) zeros(size(paper3_DV_lag)) m_t_agg(:,4) paper3_DV_lag zeros(size(paper3_DV_lag)) m_t_agg(:,4).*paper3_DV_lag];
    y_pred = SimuX*modelresults;

    for t = 1:10
        APOP_ready_vec_rec = []
        lag_y_rec = []
        for i = 1:length(m_list_rec)
            lag_y_temp = y_pred(m_start_rec(i):m_end_rec(i)-t);
            APOP_temp = APOP_ready_vec(m_start_rec(i)+t:m_end_rec(i),:);
            APOP_ready_vec_rec = [APOP_ready_vec_rec; APOP_temp];
            lag_y_rec = [lag_y_rec; lag_y_temp];
        end
        SimuX_rec = [ones(size(lag_y_rec,1),1) zeros(size(lag_y_rec,1),1) APOP_ready_vec_rec lag_y_rec zeros(size(lag_y_rec,1),1) APOP_ready_vec_rec.*lag_y_rec];
        y_rec = SimuX_rec*modelresults;
    end
    simu_paper3{s} = y_rec;
end

% save('simu_paper3.mat','simu_paper3');

load('simu_paper3.mat')

simu_paper3_temp = cell2mat(simu_paper3);
simu_paper3_mat = reshape(simu_paper3_temp,17931,300);

ww = mean(simu_paper3_mat,1);
mean(ww)
std(ww)

mean(paper3_DV)

% w = mean(simu_paper3_mat,1);

(mean(paper3_DV)-mean(ww))/mean(paper3_DV)

std(paper3_DV)
