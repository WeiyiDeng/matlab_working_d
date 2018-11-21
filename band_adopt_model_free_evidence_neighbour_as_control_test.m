clc
clear

%% source(member) & neighbour

adoptions = csvread('bandadoptions_lenient_adopt.csv',1,0);
% timesplit = csvread('tsplit3.csv');
friendlist = csvread('new_friendlist_8088.csv',1,0);
% innov = csvread('EAi3_lenient.csv');
adoptions_neighbour = csvread('bandadoptions_neighbour.csv',1,0);
neighbourlist = csvread('neighbourlist_6585.csv',1,0);        % 6585 member & neighbour pairs

adoptions_full = [adoptions(:,1:3);adoptions_neighbour];      % add neighbour adoptions to below

friendlist = [friendlist(end,:);friendlist];
friendlist(end,:) = [];

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

row_start_member2 = zeros(size(friendlist,1),1);
row_end_member2 = zeros(size(friendlist,1),1);
row_start_friend = zeros(size(friendlist,1),1);
row_end_friend = zeros(size(friendlist,1),1);
for i = 1:size(friendlist,1)
    source_index = find(user_id_store==friendlist(i,1));
    recepient_index = find(user_id_store==friendlist(i,2));
    if length(source_index)>1 || length(recepient_index)>1     % remove replicated records in adoptions_full
        source_index = source_index(1);
        recepient_index = recepient_index(1);
    end
    row_start_member2(i) = row_start(source_index);
    row_end_member2(i) = row_end(source_index);
    row_start_friend(i) = row_start(recepient_index);
    row_end_friend(i) = row_end(recepient_index);
end

%%
bands_shared = zeros(size(neighbourlist,1),1);
bands_shared_temp_neighbour = zeros(size(neighbourlist,1),1);
% bands_adopted_within_nWeeks = zeros(size(neighbourlist,1),1);
count_member_adopts = zeros(size(neighbourlist,1),1);
count_neighbour_adopts = zeros(size(neighbourlist,1),1);
band_specific_count_shared_neighbour = zeros(6046,1);
band_specific_count_shared_neighbour_sequential52 = zeros(6046,1);
band_specific_count_shared_neighbour_sequential12 = zeros(6046,1);
band_specific_count_shared_neighbour_sequential4 = zeros(6046,1);
% length_both_active_period_weeks = zeros(size(neighbourlist,1),1);
for i = 1:size(neighbourlist,1)
    source_adopts_bandIDs = adoptions_full(row_start_member(i):row_end_member(i),2);
    recepient_adopts_bandIDs = adoptions_full(row_start_neighbour(i):row_end_neighbour(i),2);
    source_adopts_weekIDs = adoptions_full(row_start_member(i):row_end_member(i),3);
    recepient_adopts_weekIDs = adoptions_full(row_start_neighbour(i):row_end_neighbour(i),3);
    for j = 1:length(source_adopts_bandIDs)
        if ismember(source_adopts_bandIDs(j),recepient_adopts_bandIDs)
            bands_shared(i) = bands_shared(i)+1;
        end
    end
    temp1 = ismember(1:6046,source_adopts_bandIDs);
    temp2 = ismember(1:6046,recepient_adopts_bandIDs);
    bands_shared_temp_neighbour(i) = sum(temp1.*temp2);
    band_specific_count_shared_neighbour = band_specific_count_shared_neighbour + (temp1.*temp2)';
    temp_sparse = sparse(1,source_adopts_bandIDs,source_adopts_weekIDs,1,6046);
    temp_sparse2 = sparse(1,recepient_adopts_bandIDs,recepient_adopts_weekIDs,1,6046);
    logical_index52 = temp_sparse.*temp_sparse2>0 & abs(temp_sparse-temp_sparse2)<52;
    logical_index12 = temp_sparse.*temp_sparse2>0 & abs(temp_sparse-temp_sparse2)<12;
    logical_index4 = temp_sparse.*temp_sparse2>0 & abs(temp_sparse-temp_sparse2)<4;
    band_specific_count_shared_neighbour_sequential52 = band_specific_count_shared_neighbour_sequential52+logical_index52';
    band_specific_count_shared_neighbour_sequential12 = band_specific_count_shared_neighbour_sequential12+logical_index12';
    band_specific_count_shared_neighbour_sequential4 = band_specific_count_shared_neighbour_sequential4+logical_index4';
    count_member_adopts(i) = length(source_adopts_bandIDs);
    count_neighbour_adopts(i) = length(recepient_adopts_bandIDs);
end

sum(band_specific_count_shared_neighbour) == sum(bands_shared_temp_neighbour)

                    % threshold of sequential behaviour
bands_shared_f = zeros(size(friendlist,1),1);
bands_shared_temp_friend = zeros(size(friendlist,1),1);
% bands_adopted_within_nWeeks = zeros(size(neighbourlist,1),1);
count_member_adopts2 = zeros(size(friendlist,1),1);
count_friend_adopts = zeros(size(friendlist,1),1);
band_specific_count_shared_friend = zeros(6046,1);
band_specific_count_shared_friend_sequential52 = zeros(6046,1);
band_specific_count_shared_friend_sequential12 = zeros(6046,1);
band_specific_count_shared_friend_sequential4 = zeros(6046,1);
% length_both_active_period_weeks = zeros(size(neighbourlist,1),1);
for i = 1:size(friendlist,1)
    source_adopts_bandIDs = adoptions_full(row_start_member2(i):row_end_member2(i),2);
    recepient_adopts_bandIDs = adoptions_full(row_start_friend(i):row_end_friend(i),2);
    source_adopts_weekIDs = adoptions_full(row_start_member2(i):row_end_member2(i),3);
    recepient_adopts_weekIDs = adoptions_full(row_start_friend(i):row_end_friend(i),3);
    for j = 1:length(source_adopts_bandIDs)
        if ismember(source_adopts_bandIDs(j),recepient_adopts_bandIDs)
            bands_shared_f(i) = bands_shared_f(i)+1;
        end
    end
    temp1 = ismember(1:6046,source_adopts_bandIDs);
    temp2 = ismember(1:6046,recepient_adopts_bandIDs);
    bands_shared_temp_friend(i) = sum(temp1.*temp2);
    band_specific_count_shared_friend = band_specific_count_shared_friend + (temp1.*temp2)';
    temp_sparse = sparse(1,source_adopts_bandIDs,source_adopts_weekIDs,1,6046);
    temp_sparse2 = sparse(1,recepient_adopts_bandIDs,recepient_adopts_weekIDs,1,6046);
    logical_index52 = temp_sparse.*temp_sparse2>0 & abs(temp_sparse-temp_sparse2)<52;
    logical_index12 = temp_sparse.*temp_sparse2>0 & abs(temp_sparse-temp_sparse2)<12;
    logical_index4 = temp_sparse.*temp_sparse2>0 & abs(temp_sparse-temp_sparse2)<4;
    band_specific_count_shared_friend_sequential52 = band_specific_count_shared_friend_sequential52+logical_index52';
    band_specific_count_shared_friend_sequential12 = band_specific_count_shared_friend_sequential12+logical_index12';
    band_specific_count_shared_friend_sequential4 = band_specific_count_shared_friend_sequential4+logical_index4';
    count_member_adopts2(i) = length(source_adopts_bandIDs);
    count_friend_adopts(i) = length(recepient_adopts_bandIDs);
end

sum(band_specific_count_shared_friend) == sum(bands_shared_temp_friend)

% normalized for # of friends & # ofneighbours
figure
threshold = 400;
keep_index = band_specific_count_shared_friend<threshold & band_specific_count_shared_neighbour<threshold;

% x = log10(1+band_specific_count_shared_friend(keep_index)./size(friendlist,1));
% y = log10(1+band_specific_count_shared_neighbour(keep_index)./size(neighbourlist,1));
% x = band_specific_count_shared_friend_sequential52(keep_index)./size(friendlist,1);
% y = band_specific_count_shared_neighbour_sequential52(keep_index)./size(neighbourlist,1);
% x = band_specific_count_shared_friend_sequential12(keep_index)./size(friendlist,1);
% y = band_specific_count_shared_neighbour_sequential12(keep_index)./size(neighbourlist,1);
x = log10(1+band_specific_count_shared_friend_sequential4(keep_index)./size(friendlist,1));
y = log10(1+band_specific_count_shared_neighbour_sequential4(keep_index)./size(neighbourlist,1));

scatter(x, y,'.')
hold on

Fit = polyfit(x,y,1);
xx = [min(x) max(y)];
yhat = polyval(Fit,xx);
plot(xx,yhat,'r--')
hold on

plot([0 0.006],[0 0.006],'color',[230, 159, 0] / 256)
hold off
xlabel('shared with friends') 
ylabel('shared with neigbours')
title('bands shared with # of friends vs # of neighbours normalized by number of friends/neighours')

%%
total_bands_N = 440000;
% total_bands_N = 4400
BP_MN = count_member_adopts.*count_neighbour_adopts/total_bands_N;
X_MN = bands_shared./BP_MN;

BP_MF = count_member_adopts2.*count_friend_adopts/total_bands_N;
X_MF = bands_shared_f./BP_MF;
% notice X_MF(2822), friendlist(2822,:),
% find(members_in_friendlist==5567)=>81
% find(members_in_neighbourlist==5567)=>78
% members_row_start_friendlist(81), members_row_end_friendlist(81)
% members_row_start_neighbourlist(78), members_row_end_neighbourlist(78)

members_in_neighbourlist = [];
members_row_start_neighbourlist = 1;
members_row_end_neighbourlist = [];
member_id_store = neighbourlist(1,1);
for i=2:size(neighbourlist,1)
    member_id = neighbourlist(i,1);
    if member_id~=member_id_store
        members_in_neighbourlist = [members_in_neighbourlist member_id_store];
        members_row_end_neighbourlist = [members_row_end_neighbourlist i-1];
        members_row_start_neighbourlist = [members_row_start_neighbourlist i];
    end
    member_id_store = member_id;
end
members_in_neighbourlist = [members_in_neighbourlist neighbourlist(end,1)];
members_row_end_neighbourlist = [members_row_end_neighbourlist size(neighbourlist,1)];

members_in_friendlist = [];
members_row_start_friendlist = 1;
members_row_end_friendlist = [];
member_id_store = friendlist(1,1);
for i=2:size(friendlist,1)
    member_id = friendlist(i,1);
    if member_id~=member_id_store
        members_in_friendlist = [members_in_friendlist member_id_store];
        members_row_end_friendlist = [members_row_end_friendlist i-1];
        members_row_start_friendlist = [members_row_start_friendlist i];
    end
    member_id_store = member_id;
end
members_in_friendlist = [members_in_friendlist friendlist(end,1)];
members_row_end_friendlist = [members_row_end_friendlist size(friendlist,1)];

MN_MF_combined = zeros(length(members_in_neighbourlist),2);
avg_shared_MN = [];
avg_shared_MF = [];
sum(ismember(members_in_neighbourlist,members_in_friendlist))     % all members_in_neighbourlist appear in members_in_friendlist
for j=1:length(members_in_neighbourlist)
    member_index = find(members_in_neighbourlist(j)==members_in_friendlist);
    member_neighbour_rows = members_row_start_neighbourlist(j):members_row_end_neighbourlist(j);
    member_friend_rows = members_row_start_friendlist(member_index):members_row_end_friendlist(member_index);
    MN_MF_combined(j,1) = mean(X_MN(member_neighbour_rows));
    MN_MF_combined(j,2) = mean(X_MF(member_friend_rows));
    avg_shared_MN = mean(X_MN(member_neighbour_rows).*BP_MN(member_neighbour_rows));
    avg_shared_MF = mean(X_MN(member_neighbour_rows).*BP_MN(member_neighbour_rows));
end

%%
NF_bands_shared = zeros(500000,1);
NF_count_bands_neighbour_adopts = zeros(500000,1);
NF_count_bands_friend_adopts = zeros(500000,1);
NF_member_index = zeros(500000,1);
avg_BP_NF = [];
avg_X_NF = [];
avg_shared_NF = [];
count_ind = 0;
for n=1:length(members_in_neighbourlist)
    member_index = find(members_in_neighbourlist(n)==members_in_friendlist);
    member_neighbour_rows = members_row_start_neighbourlist(n):members_row_end_neighbourlist(n);
    member_friend_rows = members_row_start_friendlist(member_index):members_row_end_friendlist(member_index);
    neighbours_ind = neighbourlist(member_neighbour_rows,2);
    friends_ind = friendlist(member_friend_rows,2);
    per_member_BP_NF = [];
    per_member_X_NF = [];
    per_member_shared_NF = [];
    for i =1:length(neighbours_ind)
        neighbour_adopts_ind = find(neighbours_ind(i)==user_id_store);
        neighbour_adopts = row_start(neighbour_adopts_ind):row_end(neighbour_adopts_ind);
        bands_neighbour_adopts = adoptions_full(neighbour_adopts,2);
        for j = 1:length(friends_ind)
            friend_adopts_ind = find(friends_ind(j)==user_id_store);
            friend_adopts = row_start(friend_adopts_ind):row_end(friend_adopts_ind);
            bands_friend_adopts = adoptions_full(friend_adopts,2);
            bands_sharedij = 0;
            for k = 1:length(bands_neighbour_adopts)
                if ismember(bands_neighbour_adopts(k),bands_friend_adopts)
                    bands_sharedij = bands_sharedij+1;
                end
            end
            count_ind = count_ind+1;
            NF_bands_shared(count_ind) = bands_sharedij;
            NF_count_bands_neighbour_adopts(count_ind) = length(bands_neighbour_adopts);
            NF_count_bands_friend_adopts(count_ind) = length(bands_friend_adopts);
            NF_member_index(count_ind) = member_index;
%             NF_bands_shared = [NF_bands_shared bands_sharedij];
%             NF_count_bands_neighbour_adopts = [NF_count_bands_neighbour_adopts length(bands_neighbour_adopts)];
%             NF_count_bands_friend_adopts = [NF_count_bands_friend_adopts length(bands_friend_adopts)];
%             NF_member_index = [NF_member_index member_index];
            BP_NF = length(bands_neighbour_adopts)*length(bands_friend_adopts)/total_bands_N;
            X_NF = bands_sharedij/BP_NF;
            per_member_BP_NF = [per_member_BP_NF BP_NF];
            per_member_X_NF = [per_member_X_NF X_NF];
            per_member_shared_NF = [per_member_shared_NF bands_sharedij];
        end
    end
    avg_BP_NF = [avg_BP_NF mean(per_member_BP_NF)];
    avg_X_NF = [avg_X_NF mean(per_member_X_NF)];
    avg_shared_NF = [avg_shared_NF mean(per_member_shared_NF)];
end

NF_bands_shared = NF_bands_shared(1:count_ind);
NF_count_bands_neighbour_adopts = NF_count_bands_neighbour_adopts(1:count_ind);
NF_count_bands_friend_adopts = NF_count_bands_friend_adopts(1:count_ind);
NF_member_index = NF_member_index(1:count_ind);

avg_X_MF = MN_MF_combined(:,2);
compareX_MF_NF = avg_X_MF-avg_X_NF';

% w: if # of band adoptions is small, but have some popular bands shared by
% both parties, then X scaled by BP will be blown up
% w: do neighbours tend to adopt more? Yes, but not by a large proportion
% (30%) see below
% w: try t test 

%% test: do neighbours tend to adopt more
% neighbour_ids = unique(neighbourlist(:,2));
% friend_ids = unique(friendlist(:,2));
% 
% count_neighbours_adopt = [];
% for i = 1:length(neighbour_ids)
%     rows =find(adoptions_full(:,1)==neighbour_ids(i));
%     count_neighbours_adopt = [count_neighbours_adopt length(rows)];
% end
% 
% count_friends_adopt = [];
% for i = 1:length(friend_ids)
%     rows =find(adoptions_full(:,1)==friend_ids(i));
%     count_friends_adopt = [count_friends_adopt length(rows)];
% end
% 
% mean(count_neighbours_adopt)   % => 134
% mean(count_friends_adopt)      % => 100
% hist(count_neighbours_adopt)
% hist(count_friends_adopt)

%% run Wilcoxon rank sum test (equivalent to a Mann-Whitney U-test)
% https://www.mathworks.com/help/stats/ranksum.html?s_tid=gn_loc_drop
avg_X_NF = avg_X_NF';
[p,h,stats] = ranksum(avg_X_MF,avg_X_NF) 

% [p,h,stats] = ranksum(avg_shared_MF,avg_shared_NF')

% visualize in histogram
a = avg_X_MF-avg_X_NF;
hist(a(a<1000))


