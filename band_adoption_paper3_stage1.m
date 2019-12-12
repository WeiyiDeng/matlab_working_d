clc
clear

date_=clock;
resultsfilename=['Results/r_Results_' num2str(date_(1)) '_' num2str(date_(2)) '_' num2str(date_(3)) '-' num2str(date_(4))  '_' num2str(date_(5)) '.txt'];
diary(resultsfilename);

mfilename('fullpath')                           % print filename of currently running script

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
load('dummy_mat.mat');

matp = matp(indx,:);

X = dummy_mat(:, 5);
y = dummy_mat(:,4);
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

clearvars N_prep_for_matp_jobs_organize combi_similarities matp

%%
% load('b_agg_dummy_pop.mat')
% beta_0 = b
beta_0 = [-7.7250    0.0686   0.01    0.01    0.01    0.001   0.0341    2.1842   0.7598]
% beta_0 = [-7.7250    0.0686   0.01    0.01    0.0341    2.1842   0.7598]
% beta_0 = [-7.7250    0.0686   0.01    0.0341    2.1842   0.7598]

clearvars b

[b, hessian, grad, standard_error, covariance_matrix, t_stat, exit_flag, output] = band_runbi_ll_paper3_stage1(X, trend_hat, pop, dummy_agg_SI, None0s_X_N, S, y, beta_0);

% save('b_paper3_stage1_noGT.mat','b') ;
% save('standard_error_paper3_stage1_noGT.mat','standard_error') ;
% save('t_stat_paper3_stage1_noGT.mat','t_stat') ;
% save('exit_flag_paper3_stage1_noGT.mat','exit_flag') ;
% save('hessian_paper3_stage1_noGT.mat','hessian') ;

% save('b_paper3_stage1_noGT_noInteract.mat','b') ;
% save('standard_error_paper3_stage1_noGT_noInteract.mat','standard_error') ;
% save('t_stat_paper3_stage1_noGT_noInteract.mat','t_stat') ;
% save('exit_flag_paper3_stage1_noGT_noInteract.mat','exit_flag') ;
% save('hessian_paper3_stage1_noGT_noInteract.mat','hessian') ;

% save('b_paper3_stage1_noGT_noInteract_noQrad.mat','b') ;
% save('standard_error_paper3_stage1_noGT_noInteract_noQrad.mat','standard_error') ;
% save('t_stat_paper3_stage1_noGT_noInteract_noQrad.mat','t_stat') ;
% save('exit_flag_paper3_stage1_noGT_noInteract_noQrad.mat','exit_flag') ;
% save('hessian_paper3_stage1_noGT_noInteract_noQrad.mat','hessian') ;

display(b)
% display(standard_error)
display(t_stat)
display(grad)
display(output)

diary off 


