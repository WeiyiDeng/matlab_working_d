clc
clear

adoptions = csvread('bandadoptions_lenient_adopt.csv',1,0);
timesplit = csvread('tsplit3.csv');
friendlist = csvread('new_friendlist_8088.csv',1,0);
innov = csvread('EAi3_lenient.csv');

temp_tsplit = zeros(size(adoptions,1),1);
for i = 1:length(temp_tsplit)
    user_id = adoptions(i,1);
    temp_tsplit(i) = timesplit(user_id,3);
end

temp_innov = zeros(size(adoptions,1),1);
for i = 1:length(temp_innov)
    user_id = adoptions(i,1);
    temp_innov(i) = innov(user_id,2);
end

sum(isnan(temp_innov))      % test

% denote 1 to the 2nd part adoptions observations
adoptions_2nd_half = adoptions(:,3)>temp_tsplit;
sum(adoptions_2nd_half)           % test

% set thresholds to differentiate low medium and high adoptions
l_innov = quantile(temp_innov,0.33)         
m_innov = quantile(temp_innov,0.66)
% l_innov = quantile(temp_innov,0.5)         
% m_innov = quantile(temp_innov,0.9)

discrete_innov = zeros(size(adoptions,1),1);
for i = 1:length(discrete_innov)
    if temp_innov(i) < l_innov
        discrete_innov(i) = 1;             % low innovativeness
    elseif temp_innov(i) > m_innov
        discrete_innov(i) = 3;             % high
    else
        discrete_innov(i) = 2;             % medium
    end
end

temp_user_id = 2;
row_start = 1;
row_end = [];
for i = 1:size(adoptions,1)
    if adoptions(i,1)==temp_user_id
    else
        row_end = [row_end i-1];
        row_start = [row_start i];
        temp_user_id = adoptions(i,1);
    end
end
row_end = [row_end i];
user_id_store = adoptions(row_start,1);
        
count_2nd_half_adoptions = zeros(size(user_id_store));
for i = 1:length(user_id_store)
    count_2nd_half_adoptions(i) = sum(adoptions_2nd_half(row_start(i):row_end(i)));
end
discrete_innov_short = discrete_innov(row_start); 
    
% hist(count_2nd_half_adoptions,50)    

% csvwrite('count_2nd_half_adoptions.csv',count_2nd_half_adoptions);
% csvwrite('discrete_innov_short.csv',discrete_innov_short);
    
%% source & recepient 
% (w: not using 2nd half observations at the moment !!)
source_innov_level = zeros(size(friendlist,1),1);
recepient_innov_level = zeros(size(friendlist,1),1);
for i = 1:size(friendlist,1)
    row_index_s = find(user_id_store==friendlist(i,1));
    row_index_r = find(user_id_store==friendlist(i,2));
    source_innov_level(i) = discrete_innov_short(row_index_s);
    recepient_innov_level(i) = discrete_innov_short(row_index_r);
end
HH_index = source_innov_level==1 & recepient_innov_level==2;
sum(HH_index)    
    
% HH_friendlist = friendlist(HH_index,:);                  % select friends and members pairs both with high innovative scores
HH_friendlist = friendlist;                      % not distinguishing between different innovativeness types (select all users for now)

HH_row_start_s = zeros(size(HH_friendlist,1),1);
HH_row_end_s = zeros(size(HH_friendlist,1),1);
HH_row_start_r = zeros(size(HH_friendlist,1),1);
HH_row_end_r = zeros(size(HH_friendlist,1),1);
for i = 1:size(HH_friendlist,1)
    source_index = find(user_id_store==HH_friendlist(i,1));
    recepient_index = find(user_id_store==HH_friendlist(i,2));
    HH_row_start_s(i) = row_start(source_index);
    HH_row_end_s(i) = row_end(source_index);
    HH_row_start_r(i) = row_start(recepient_index);
    HH_row_end_r(i) = row_end(recepient_index);
end

n = 4;            % set within n number of weeks
bands_shared = zeros(size(HH_friendlist,1),1);
bands_adopted_within_nWeeks = zeros(size(HH_friendlist,1),1);
count_source_adopts = zeros(size(HH_friendlist,1),1);
count_recepient_adopts = zeros(size(HH_friendlist,1),1);
length_both_active_period_weeks = zeros(size(HH_friendlist,1),1);
for i = 1:size(HH_friendlist,1)
    source_adopts_bandIDs = adoptions(HH_row_start_s(i):HH_row_end_s(i),2);
    recepient_adopts_bandIDs = adoptions(HH_row_start_r(i):HH_row_end_r(i),2);
    source_adopts_weekIDs = adoptions(HH_row_start_s(i):HH_row_end_s(i),3);
    recepient_adopts_weekIDs = adoptions(HH_row_start_r(i):HH_row_end_r(i),3);
    source_start_platform_active = min(source_adopts_weekIDs);               % find active period on platform of source and recepient
    source_end_platform_active = max(source_adopts_weekIDs);
    recepient_start_platform_active = min(recepient_adopts_weekIDs);
    recepient_end_platform_active = max(recepient_adopts_weekIDs);
    lower_bound_active = max(source_start_platform_active,recepient_start_platform_active);
    upper_bound_active = min(source_end_platform_active,recepient_end_platform_active);
    if lower_bound_active > upper_bound_active                               
        continue                             % for source & recepient with no overlap active period, skip and move to the next iteration                 
    end
    same_period_source_adopts_bandIDs = source_adopts_bandIDs(source_adopts_weekIDs>lower_bound_active &...
        source_adopts_weekIDs<upper_bound_active);
    same_period_recepient_adopts_bandIDs = recepient_adopts_bandIDs(recepient_adopts_weekIDs>lower_bound_active &...
        recepient_adopts_weekIDs<upper_bound_active);
    same_period_source_adopts_weekIDs = source_adopts_weekIDs(source_adopts_weekIDs>lower_bound_active &...
        source_adopts_weekIDs<upper_bound_active);
    same_period_recepient_adopts_weekIDs = recepient_adopts_weekIDs(recepient_adopts_weekIDs>lower_bound_active &...
        recepient_adopts_weekIDs<upper_bound_active);
    count_source_adopts(i) = numel(same_period_source_adopts_bandIDs);           % # of adoptions by source within active period on platform by both
    count_recepient_adopts(i) = numel(same_period_recepient_adopts_bandIDs);     % # of adoptions by recepient within active period on platform by both
    length_both_active_period_weeks(i) = upper_bound_active-lower_bound_active;
    for j = 1:length(same_period_source_adopts_bandIDs)
        if ismember(same_period_source_adopts_bandIDs(j),same_period_recepient_adopts_bandIDs)
            bands_shared(i) = bands_shared(i)+1;                                  % shared band adoptions within active period on platform by both
            recepient_band_index = find(same_period_recepient_adopts_bandIDs==same_period_source_adopts_bandIDs(j));
            if same_period_recepient_adopts_weekIDs(recepient_band_index)-same_period_source_adopts_weekIDs(j)< n && ...  % within n number of weeks
                    same_period_recepient_adopts_weekIDs(recepient_band_index)-same_period_source_adopts_weekIDs(j)>0
                bands_adopted_within_nWeeks(i) = bands_adopted_within_nWeeks(i)+1;
            else
            end
        else
        end
    end
end

count_sum_source_recepient_adopts = count_source_adopts+count_recepient_adopts;
bands_shared = bands_shared(count_sum_source_recepient_adopts~=0);
count_sum_source_recepient_adopts = count_sum_source_recepient_adopts(count_sum_source_recepient_adopts~=0);
homophily_index = bands_shared./count_sum_source_recepient_adopts;
% disp(mean(homophily_index))
disp(['homophily_index    ' num2str(mean(homophily_index))])

bands_adopted_within_nWeeks = bands_adopted_within_nWeeks(count_sum_source_recepient_adopts~=0);
simultaneous_behavior = bands_adopted_within_nWeeks./count_sum_source_recepient_adopts;
% disp(mean(simultaneous_behavior))
% 
% length_both_active_period_weeks = length_both_active_period_weeks(count_sum_source_recepient_adopts~=0);
% length_both_active_period_weeks = length_both_active_period_weeks(length_both_active_period_weeks~=0);
% simultaneous_behavior = simultaneous_behavior(length_both_active_period_weeks~=0);
% homophily_index = homophily_index(length_both_active_period_weeks~=0);
% social_impact = simultaneous_behavior - homophily_index.*n./length_both_active_period_weeks;
% disp(mean(social_impact))

disp(['simultaneous_behavior    ' num2str(mean(simultaneous_behavior))])
% disp(['social_impact    ' num2str(mean(social_impact))])

