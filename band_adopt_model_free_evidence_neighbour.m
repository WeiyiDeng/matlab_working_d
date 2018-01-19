clc
clear

%% source(member) & neighbour

adoptions = csvread('bandadoptions_lenient_adopt.csv',1,0);
% timesplit = csvread('tsplit3.csv');
friendlist = csvread('new_friendlist_8088.csv',1,0);
% innov = csvread('EAi3_lenient.csv');
adoptions_neighbour = csvread('bandadoptions_neighbour.csv',1,0);
neighbourlist = csvread('neighbourlist_6585.csv',1,0);        % 6585 member & neighbour pairs

adoptions_full = [adoptions(:,1:3);adoptions_neighbour];

temp_user_id = 2;
row_start = 1;
row_end = [];
for i = 1:size(adoptions_full,1)
    if adoptions_full(i,1)==temp_user_id
    else
        row_end = [row_end i-1];
        row_start = [row_start i];
        temp_user_id = adoptions_full(i,1);
    end
end
row_end = [row_end i];
user_id_store = adoptions_full(row_start,1);

sum(1-ismember(neighbourlist(:,2),user_id_store))   % 170 neighbours not found in adoptions_neighbour (no listen records)
neighbourlist(find(1-ismember(neighbourlist(:,2),user_id_store)),:) = []; % 6415 member & neighbour pairs left
sum(1-ismember(neighbourlist(:,1),user_id_store))   % 21 members not found in adoptions_neighbour (no listen records)
neighbourlist(find(1-ismember(neighbourlist(:,1),user_id_store)),:) = []; % 6394 member & neighbour pairs left

row_start_member = zeros(size(neighbourlist,1),1);
row_end_member = zeros(size(neighbourlist,1),1);
row_start_neighbour = zeros(size(neighbourlist,1),1);
row_end_neighbour = zeros(size(neighbourlist,1),1);
for i = 1:size(neighbourlist,1)
    source_index = find(user_id_store==neighbourlist(i,1));
    recepient_index = find(user_id_store==neighbourlist(i,2));
    if length(source_index)>1 || length(recepient_index)>1     % remove replicated records in adoptions_full
        source_index = source_index(1);
        recepient_index = recepient_index(1);
    end
    row_start_member(i) = row_start(source_index);
    row_end_member(i) = row_end(source_index);
    row_start_neighbour(i) = row_start(recepient_index);
    row_end_neighbour(i) = row_end(recepient_index);
end

n = 4;            % set within n number of weeks
bands_shared = zeros(size(neighbourlist,1),1);
bands_adopted_within_nWeeks = zeros(size(neighbourlist,1),1);
count_source_adopts = zeros(size(neighbourlist,1),1);
count_recepient_adopts = zeros(size(neighbourlist,1),1);
length_both_active_period_weeks = zeros(size(neighbourlist,1),1);
for i = 1:size(neighbourlist,1)
    source_adopts_bandIDs = adoptions_full(row_start_member(i):row_end_member(i),2);
    recepient_adopts_bandIDs = adoptions_full(row_start_neighbour(i):row_end_neighbour(i),2);
    source_adopts_weekIDs = adoptions_full(row_start_member(i):row_end_member(i),3);
    recepient_adopts_weekIDs = adoptions_full(row_start_neighbour(i):row_end_neighbour(i),3);
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