clc
clear

% load('matpstd2_all_member.mat');
load('matp_lenient_adopt_newBP_clean_std.mat')
% load('matp_lenient_adopt.mat')

% predict_trend = csvread('predict_trend.csv',1,0);       % reads from row 2 of the doc
% predict_trend = csvread('predict_trend_log.csv',1,0);       % reads from row 2 of the doc
predict_trend = csvread('predict_trend_log4061_lenient.csv',1,0);       % reads from row 2 of the doc

% trend_hat = [];
% for i = 1:54921
% % for i = 1:size(matp,1)
%     band_num = matp(i,3);
%     week_num = matp(i,4);
%     trend_row = find(predict_trend(:,1) == band_num & predict_trend(:,2) == week_num);
%     trend_obs = predict_trend(trend_row, 3);
%     trend_hat = [trend_hat; trend_obs];
% end

check_continuous = [];
store_band_matp = [];
band_start_row_matp = [];
band_end_row_matp = [];
last_bandnum = matp(1,3);
week_start = matp(1,4);
k = 0;
for i = 1:size(matp,1)
% for i = 1:54921
    if matp(i,3)==last_bandnum
        k = k+1;
    else
        week_end = matp(i-1,4);
        week_interval = week_end - week_start + 1;
        check_continuous = [check_continuous week_interval==k];
        store_band_matp = [store_band_matp matp(i-1,3)];
        band_start_row_matp = [band_start_row_matp i-k];
        band_end_row_matp = [band_end_row_matp i-1];
        week_start = matp(i,4);
        k = 1;
        last_bandnum = matp(i,3);
    end
end
store_band_matp = [store_band_matp matp(i,3)];
band_start_row_matp = [band_start_row_matp i-k+1];
band_end_row_matp = [band_end_row_matp i];   

% distinct_bands_in_matp = unique(store_band_matp);
distinct_bands_in_matp = unique(matp(:,3));           % 4061 distinct bands in new data matrix matpstd2_all_member.mat
csvwrite('newbandsNo4061.csv',distinct_bands_in_matp);

store_band_trend = [];
last_band_ind = predict_trend(1,1);
for j = 1: size(predict_trend,1)
    if predict_trend(j,1) == last_band_ind
    else
        store_band_trend = [store_band_trend predict_trend(j-1,1)];
        last_band_ind = predict_trend(j,1);
    end
end
store_band_trend = [store_band_trend predict_trend(size(predict_trend,1),1)];
% X = predict_trend((3*423+1):4*423,:);        

trend_hat = zeros(size(matp,1),1);
for i = 1:length(store_band_matp)
% for i = 1:size(matp,1)
    band_num_matp = store_band_matp(i);
    row_start = band_start_row_matp(i);
    row_end = band_end_row_matp(i);
    trend_ind = find(store_band_trend ==band_num_matp);          % need to get google trends of more bands > 2588 !!
    trend_rows = predict_trend(((trend_ind-1)*423+1):trend_ind*423,:);
    week_nums_matp = matp(row_start:row_end, 4);
    trend_band = zeros(length(week_nums_matp),1);
    for j = 1:length(week_nums_matp)
        trend_band(j) = trend_rows(week_nums_matp(j),3);
    end
    trend_hat(row_start:row_end) = trend_band;
end

% load('matp2.mat')                          % unstandardized baseline prob

baseline_prob_store = matp(:,6);             % temp store the baseline probability in the old model

% matp(:,6) = trend_hat;

clearvars matp

%% get rid of bands without search trend data
% load('innov_contin_std2.mat');
% load('explor_contin_std2.mat');
% innov_contin = innov_contin(find(matp(:,6))>0,:);
% explor_contin = explor_contin(find(matp(:,6))>0,:);
% save('innov_contin_trend_rm.mat','innov_contin','-v7.3');
% save('explor_contin_trend_rm.mat','explor_contin','-v7.3');

% load('innov_contin_std2_all_member.mat');
load('innov_contin2_lenient_std.mat');
innov_contin = innov_contin(find(trend_hat>0),:);
save('innov_contin_trend_rm_all_member_lenient.mat','innov_contin','-v7.3');
clearvars matp

% load('explor_contin_std2_all_member.mat');
load('explor_contin2_lenient_std.mat');
explor_contin = explor_contin(find(trend_hat>0),:);
save('explor_contin_trend_rm_all_member_lenient.mat','explor_contin','-v7.3');
clearvars matp

index_keep = find(trend_hat>0);          % store indices of obs not to remove

baseline_prob_store = baseline_prob_store(find(trend_hat>0),:);
save('newbp_store_rm_all_member_lenient.mat','baseline_prob_store','-v7.3');

load('matp_lenient_adopt_newBP_clean_std.mat')
matp(:,6) = trend_hat;
matp = matp(find(trend_hat>0),:);

% %% replace missing research trend obs with mean
% mean(matp(:,6))
% sample_mean = mean(matp(find(matp(:,6))>0,6));
% ind_missing = matp(:,6)==0;
% matp(ind_missing,6) = sample_mean;
% mean(matp(:,6))

%% standardize predicted trend
matp(:,6) = (matp(:,6)-mean(matp(:,6)))./std(matp(:,6));
% mean(matp(:,6))

%% save new matp matrix
% save('matp_trend_exp.mat','matp','-v7.3');             % take exponent back to log(y_hat) as predicted trend
% save('matp_trend.mat','matp','-v7.3');
% save('matp_trend_rm.mat','matp','-v7.3');
% save('matp_trend_strict_rm.mat','matp','-v7.3');
save('matp_trend_rm_all_member_lenient.mat','matp','-v7.3');

% corr(baseline_prob_store,matp(:,6))                    % check corr between predicted trend and baseline prob
% scatter(baseline_prob_store,matp(:,6))

% %% prep matrix to use for band_adoption_run_trend.m for running the model with both baseline probability and trend data (30 variables in total)
% % load('matp_lenient_adopt.mat')
% % load('matp_lenient_adopt_newBP_std.mat')
% % load('matp_lenient_adopt_newBP_clean_std.mat')
% % load('matp_strict_adopt_newBP_clean_std.mat')
% load('matp_former_adopt_newBP_clean_std.mat')
% matp = matp(index_keep,:);
% % save('matp_bp_lenient_rm.mat','matp','-v7.3');
% % save('matp_newbp_lenient_rm.mat','matp','-v7.3');
% % save('matp_newbp_clean_lenient_rm.mat','matp','-v7.3');
% % save('matp_newbp_clean_strict_rm.mat','matp','-v7.3');
% save('matp_newbp_clean_former_rm.mat','matp','-v7.3');


