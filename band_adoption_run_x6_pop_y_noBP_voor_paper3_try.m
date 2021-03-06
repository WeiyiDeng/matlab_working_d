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

% [b, hessian, grad, standard_error, covariance_matrix, t_stat, exit_flag, output] = band_runbi_ll_x6_pop_y_noBP(X, trend_hat, pop, dummy_agg_SI, None0s_X_N, S, y, beta_0);
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

m_t_agg = csvread('prep_paper3_mat.csv');
m_t_usage = csvread('dt7_member_ID_sum_listen.csv',1,0);

m_t_usage(:,2) = m_t_usage(:,2)-104;
paper3_DV = zeros(size(m_t_agg,1),1);
for i = 1:size(m_t_agg,1)
    ind = find(m_t_usage(:,1)==m_t_agg(i,1) & m_t_usage(:,2)==m_t_agg(i,2));
    if length(ind)>1
        disp(i)
    elseif isempty(ind)
        paper3_DV(i) = 0;
    else
        paper3_DV(i) = m_t_usage(ind,3);
    end
end

sum(paper3_DV>0)
% corr(m_t_agg,paper3_DV)
inds = paper3_DV>0;
% corr(m_t_agg(inds,:),paper3_DV(inds))
corr(m_t_agg(:,3:7),paper3_DV)
corr(m_t_agg(inds,3:7),paper3_DV(inds))
corr(m_t_agg(:,3:7))
corr(m_t_agg(:,[3 5]),paper3_DV)

mean(paper3_DV(inds))
hist(paper3_DV,50)
sum(paper3_DV==0)
length(paper3_DV)

X = [ones(size(m_t_agg,1),1) m_t_agg(:,3:7)];
y = paper3_DV;
b = inv(X'*X)*X'*y
Var = 1/(size(X,1)-size(X,2)+1)*sum((y-X*b).^2)*inv(X'*X)
t = b./sqrt(diag(Var))

% csvwrite('paper3_X.csv',X);
% csvwrite('paper3_y.csv',y);

X = [ones(size(m_t_agg,1),1) m_t_agg(:,[3 5])];
y = paper3_DV;
b = inv(X'*X)*X'*y
Var = 1/(size(X,1)-size(X,2)+1)*sum((y-X*b).^2)*inv(X'*X)
t = b./sqrt(diag(Var))

mt_merge = [m_t_agg paper3_DV];

collect_dummy_mat = dummy_mat(find(dummy_mat(:,4)),:);
sparse_dummy = sparse(collect_dummy_mat(:,1),collect_dummy_mat(:,3),1,8320,423);
paper1_DV = zeros(size(mt_merge,1),1);
for i = 1:size(mt_merge,1)
    paper1_DV(i) = sparse_dummy(mt_merge(i,1),mt_merge(i,2));
end
paper1_DV(paper1_DV>0)=1;    
    
member = mt_merge(1,1);
m_end = [];
m_start = 1;
m_list = [];
for i = 2:size(mt_merge,1)
    prev_m = member;
    member = mt_merge(i,1);
    if prev_m ~=member
        m_end = [m_end; i-1];
        m_start = [m_start; i];
        m_list = [m_list; prev_m];
    else
    end
end
m_list = [m_list; member];    
m_end = [m_end; i];    
    
sum(m_end - m_start >2)

member_fix = prep_paper3_fix(1,1);
m_end_fix = [];
m_start_fix = 1;
m_list_fix = [];
for i = 2:size(prep_paper3_fix,1)
    prev_m = member_fix;
    member_fix = prep_paper3_fix(i,1);
    if prev_m ~=member_fix
        m_end_fix = [m_end_fix; i-1];
        m_start_fix = [m_start_fix; i];
        m_list_fix = [m_list_fix; prev_m];
    else
    end
end
m_list_fix = [m_list_fix; member_fix];    
m_end_fix = [m_end_fix; i];    
    
sum(m_end_fix - m_start_fix >2)

%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prep_paper3_fix_cut = [];
% prep_paper3_fix_paper1_DV_cut = [];
% for i = 1:length(m_list_fix)
%     prep_paper3_fix_paper1_DV_cut = [prep_paper3_fix_paper1_DV_cut; prep_paper3_fix(m_start_fix(i)+1:m_end_fix(i),9)];
%     prep_paper3_fix_cut = [prep_paper3_fix_cut; prep_paper3_fix(m_start_fix(i):m_end_fix(i)-1,1:8)];
% end
% 
% prep_paper3_fix_move = [prep_paper3_fix_cut  prep_paper3_fix_paper1_DV_cut];
% % csvwrite('prep_paper3_fix_move.csv',prep_paper3_fix_move)

%% try
prep_paper3_fix_cut = [];
prep_paper3_fix_paper1_DV_cut = [];
prep_paper3_fix_paper1_DV_cut2 = [];
prep_paper3_fix_paper1_DV_cut3 = [];
prep_paper3_fix_paper1_DV_cut4 = [];
prep_paper3_fix_paper1_DV_cut5 = [];
prep_paper3_fix_paper1_DV_cut6 = [];
prep_paper3_fix_paper1_DV_cut7 = [];
prep_paper3_fix_paper1_DV_cut8 = [];
prep_paper3_fix_paper1_DV_cut9 = [];
prep_paper3_fix_paper1_DV_cut10 = [];
for i = 1:length(m_list_fix)
    prep_paper3_fix_paper1_DV_cut = [prep_paper3_fix_paper1_DV_cut; prep_paper3_fix(m_start_fix(i)+1:m_end_fix(i)-9,9)];
    prep_paper3_fix_paper1_DV_cut2 = [prep_paper3_fix_paper1_DV_cut2; prep_paper3_fix(m_start_fix(i)+2:m_end_fix(i)-8,9)];
    prep_paper3_fix_paper1_DV_cut3 = [prep_paper3_fix_paper1_DV_cut3; prep_paper3_fix(m_start_fix(i)+3:m_end_fix(i)-7,9)];
    prep_paper3_fix_paper1_DV_cut4 = [prep_paper3_fix_paper1_DV_cut4; prep_paper3_fix(m_start_fix(i)+4:m_end_fix(i)-6,9)];
    prep_paper3_fix_paper1_DV_cut5 = [prep_paper3_fix_paper1_DV_cut5; prep_paper3_fix(m_start_fix(i)+5:m_end_fix(i)-5,9)];
    prep_paper3_fix_paper1_DV_cut6 = [prep_paper3_fix_paper1_DV_cut6; prep_paper3_fix(m_start_fix(i)+6:m_end_fix(i)-4,9)];
    prep_paper3_fix_paper1_DV_cut7 = [prep_paper3_fix_paper1_DV_cut7; prep_paper3_fix(m_start_fix(i)+7:m_end_fix(i)-3,9)];
    prep_paper3_fix_paper1_DV_cut8 = [prep_paper3_fix_paper1_DV_cut8; prep_paper3_fix(m_start_fix(i)+8:m_end_fix(i)-2,9)];
    prep_paper3_fix_paper1_DV_cut9 = [prep_paper3_fix_paper1_DV_cut9; prep_paper3_fix(m_start_fix(i)+9:m_end_fix(i)-1,9)];
    prep_paper3_fix_paper1_DV_cut10 = [prep_paper3_fix_paper1_DV_cut10; prep_paper3_fix(m_start_fix(i)+10:m_end_fix(i),9)];
    prep_paper3_fix_cut = [prep_paper3_fix_cut; prep_paper3_fix(m_start_fix(i):m_end_fix(i)-10,1:8)];
end

prep_paper3_fix_move_agg10 = [(1:length(prep_paper3_fix_paper1_DV_cut))'  prep_paper3_fix_cut];
prep_paper3_fix_paper1_DV_cut10 = [(1:length(prep_paper3_fix_paper1_DV_cut))'  prep_paper3_fix_paper1_DV_cut  prep_paper3_fix_paper1_DV_cut2...
    prep_paper3_fix_paper1_DV_cut3  prep_paper3_fix_paper1_DV_cut4  prep_paper3_fix_paper1_DV_cut5...  
    prep_paper3_fix_paper1_DV_cut6  prep_paper3_fix_paper1_DV_cut7...
    prep_paper3_fix_paper1_DV_cut8  prep_paper3_fix_paper1_DV_cut9  prep_paper3_fix_paper1_DV_cut10];

% clearvars -EXCEPT prep_paper3_fix_move_agg10 prep_paper3_fix_paper1_DV_cut10
% 
% csvwrite('prep_paper3_fix_move_agg10.csv',prep_paper3_fix_move_agg10)
% csvwrite('prep_paper3_fix_paper1_DV_cut10.csv',prep_paper3_fix_paper1_DV_cut10)
% 
% csvwrite('prep_paper3_fix_move_agg10_test.csv',prep_paper3_fix_move_agg10(1:100,:))
% csvwrite('prep_paper3_fix_paper1_DV_cut10_test.csv',prep_paper3_fix_paper1_DV_cut10(1:100,:))

%%
% SI_lag_2 = [];
% pop_lag_2 = [];
% DV_lag_1 = [];
% paper1_DV_lag_1 = [];
% for i = 1:length(m_list)
%     SI_lag_2 = [SI_lag_2; 0; 0; mt_merge(m_start(i):(m_end(i)-2),3)];
%     pop_lag_2 = [pop_lag_2; 0; 0; mt_merge(m_start(i):(m_end(i)-2),5)];
%     DV_lag_1 = [DV_lag_1; 0; mt_merge(m_start(i):(m_end(i)-1),8)];
%     paper1_DV_lag_1 = [paper1_DV_lag_1; 0; paper1_DV(m_start(i):(m_end(i)-1))];
% end
    
SI_lag_2cut = [];
pop_lag_2cut = [];
DV_lag_1cut = [];
DV_cut = [];
paper1_DV_lag1cut = [];
m_cut = [];
t_cut = [];
for i = 1:length(m_list)
    SI_lag_2cut = [SI_lag_2cut; mt_merge(m_start(i):(m_end(i)-2),3)];
    pop_lag_2cut = [pop_lag_2cut; mt_merge(m_start(i):(m_end(i)-2),5)];
    DV_lag_1cut = [DV_lag_1cut; mt_merge((m_start(i)+1):(m_end(i)-1),8)];
    DV_cut = [DV_cut; mt_merge((m_start(i)+2):m_end(i),8)];
    m_cut = [m_cut; mt_merge((m_start(i)+2):m_end(i),1)];
    t_cut = [t_cut; mt_merge((m_start(i)+2):m_end(i),2)];
    paper1_DV_lag1cut = [paper1_DV_lag1cut; paper1_DV(m_start(i)+1:m_end(i)-1)];
end

paper3_prep = [SI_lag_2cut pop_lag_2cut DV_lag_1cut DV_cut paper1_DV_lag1cut];
paper3_met_lags = [SI_lag_2cut.*paper1_DV_lag1cut pop_lag_2cut.*paper1_DV_lag1cut DV_lag_1cut DV_cut];
paper3_lags_check = [m_cut  t_cut  SI_lag_2cut  pop_lag_2cut   DV_lag_1cut   DV_cut  paper1_DV_lag1cut];
% paper3_met_lags = [SI_lag_2cut pop_lag_2cut DV_lag_1cut DV_cut];

% csvwrite('paper3_met_lags.csv',paper3_met_lags);
paper3_met_lags = csvread('paper3_met_lags.csv');

corr(paper3_met_lags)

X = [ones(size(paper3_met_lags,1),1) paper3_met_lags(:,1:3)];
y = paper3_met_lags(:,4);
b = inv(X'*X)*X'*y
Var = 1/(size(X,1)-size(X,2)+1)*sum((y-X*b).^2)*inv(X'*X)
t = b./sqrt(diag(Var))

R_squared = corr(paper3_met_lags(:,end), X*b)^2

% X = paper3_met_lags(:,1:3);
% y = paper3_met_lags(:,4);
% mdl = fitlm(X,y)

%%
mt_fix_agg = csvread('prep_paper3_fix_mat.csv');
% mt_fix_mv_agg = csvread('prep_paper3_fix_move_mat.csv');

max(mt_fix_agg(:,2))
% max(mt_fix_mv_agg(:,2))

paper3_DV_fix = zeros(size(mt_fix_agg,1),1);
for i = 1:size(mt_fix_agg,1)
    ind = find(m_t_usage(:,1)==mt_fix_agg(i,1) & m_t_usage(:,2)==mt_fix_agg(i,2));
    if length(ind)>1
        disp(i)
    elseif isempty(ind)
        paper3_DV_fix(i) = 0;
    else
        paper3_DV_fix(i) = m_t_usage(ind,3);
    end
end

mt_merge_fix = [mt_fix_agg paper3_DV_fix];

% paper3_DV_mv_fix = zeros(size(mt_fix_mv_agg,1),1);
% for i = 1:size(mt_fix_mv_agg,1)
%     ind = find(m_t_usage(:,1)==mt_fix_mv_agg(i,1) & m_t_usage(:,2)==mt_fix_mv_agg(i,2));
%     if length(ind)>1
%         disp(i)
%     elseif isempty(ind)
%         paper3_DV_mv_fix(i) = 0;
%     else
%         paper3_DV_mv_fix(i) = m_t_usage(ind,3);
%     end
% end
% 
% mt_merge_mv_fix = [mt_fix_mv_agg paper3_DV_mv_fix];
% mt_merge_fix = mt_merge_mv_fix;

member_fix = mt_merge_fix(1,1);
m_end_fix = [];
m_start_fix = 1;
m_list_fix = [];
for i = 2:size(mt_merge_fix,1)
    prev_m = member_fix;
    member_fix = mt_merge_fix(i,1);
    if prev_m ~=member_fix
        m_end_fix = [m_end_fix; i-1];
        m_start_fix = [m_start_fix; i];
        m_list_fix = [m_list_fix; prev_m];
    else
    end
end
m_list_fix = [m_list_fix; member_fix];    
m_end_fix = [m_end_fix; i];    
    
sum(m_end_fix - m_start_fix >2)

SI_lag_2cut_fix = [];
pop_lag_2cut_fix = [];
DV_lag_1cut_fix = [];
DV_cut_fix = [];
m_cut_fix = [];
t_cut_fix = [];
for i = 1:length(m_list_fix)
    SI_lag_2cut_fix = [SI_lag_2cut_fix; mt_merge_fix(m_start_fix(i):(m_end_fix(i)-2),3)];
    pop_lag_2cut_fix = [pop_lag_2cut_fix; mt_merge_fix(m_start_fix(i):(m_end_fix(i)-2),5)];
    DV_lag_1cut_fix = [DV_lag_1cut_fix; mt_merge_fix((m_start_fix(i)+1):(m_end_fix(i)-1),8)];
    DV_cut_fix = [DV_cut_fix; mt_merge_fix((m_start_fix(i)+2):m_end_fix(i),8)];
    m_cut_fix = [m_cut_fix;mt_merge_fix((m_start_fix(i)+2):m_end_fix(i),1)];
    t_cut_fix = [t_cut_fix; mt_merge_fix((m_start_fix(i)+2):m_end_fix(i),2)];
end

paper3_prep_fix = [SI_lag_2cut_fix pop_lag_2cut_fix DV_lag_1cut_fix DV_cut_fix];
paper3_met_lags_fix = [SI_lag_2cut_fix pop_lag_2cut_fix DV_lag_1cut_fix DV_cut_fix];
paper3_lags_check_fix = [m_cut_fix  t_cut_fix  SI_lag_2cut_fix  pop_lag_2cut_fix   DV_lag_1cut_fix   DV_cut_fix];

% csvwrite('paper3_met_lags_fix.csv',paper3_met_lags_fix);
% csvwrite('paper3_met_lags_mv_fix.csv',paper3_met_lags_fix);

paper3_met_lags_fix = csvread('paper3_met_lags_fix.csv');
% paper3_met_lags_fix = csvread('paper3_met_lags_mv_fix.csv');

X = [ones(size(paper3_met_lags_fix,1),1) paper3_met_lags_fix(:,1:3)];
y = paper3_met_lags_fix(:,4);

b = inv(X'*X)*X'*y
Var = 1/(size(X,1)-size(X,2)+1)*sum((y-X*b).^2)*inv(X'*X)
t = b./sqrt(diag(Var))

R_squared = corr(y, X*b)^2

% check_mv_fix = [m_cut_fix t_cut_fix paper3_met_lags_fix];

%%
mt_fix_mv_agg2 = csvread('prep_paper3_fix_move_mat_2.csv');
paper1_NA = mt_fix_mv_agg2(:,8);
paper1_NA = paper1_NA';

paper1_NA_store = [];
for i = 1:length(m_list_fix)
    paper1_NA_store = [paper1_NA_store; paper1_NA((m_start_fix(i)+1):(m_end_fix(i)-1))];
end

X = [ones(size(paper3_met_lags_fix,1),1) paper3_met_lags_fix(:,1:3) paper1_NA_store];
y = paper3_met_lags_fix(:,4);

b = inv(X'*X)*X'*y
Var = 1/(size(X,1)-size(X,2)+1)*sum((y-X*b).^2)*inv(X'*X)
t = b./sqrt(diag(Var))

R_squared = corr(y, X*b)^2

