clc
clear

adoptions = csvread('bandadoptions_lenient_adopt.csv',1,0);
adoptions = adoptions(:,1:3);

length(unique(adoptions(:,1)))
length(unique(adoptions(:,2)))
length(unique(adoptions(:,3)))

adoptions_neighbour = csvread('bandadoptions_neighbour.csv',1,0);
adoptions_neighbour = adoptions_neighbour(adoptions_neighbour(:,1)>8320,:);       % remove duplicated user_ids in adoptions (neighbour ids appeared in member/friend ids)

length(unique(adoptions_neighbour(:,1)))
length(unique(adoptions_neighbour(:,2)))
length(unique(adoptions_neighbour(:,3)))

adoptions_full = [adoptions(:,1:3);adoptions_neighbour];

length(unique(adoptions_full(:,1)))
length(unique(adoptions_full(:,2)))
length(unique(adoptions_full(:,3)))

adopt_ones = ones(size(adoptions_full, 1),1);

% I = 13727;
I = max(adoptions_full(:,1));
J = 6046;
T = 423;
band_adoption_full = sparse(adoptions_full(:,1),adoptions_full(:,2),adopt_ones,I,J);

frequency_adopted_users = sum(band_adoption_full,1);
IDF = 1-log(frequency_adopted_users./I);
IDF_weighting_matrix = diag(IDF.^2);

%%
friendlist = csvread('new_friendlist_8088.csv',1,0);
neighbourlist = csvread('neighbourlist_6585.csv',1,0);

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

sim_mf = zeros(size(friendlist,1),1);
sim_mf_random = zeros(size(friendlist,1),1);
for i = 1:size(friendlist,1)
    member_id = friendlist(i,1);
    friend_id = friendlist(i,2);
    temp_nume = band_adoption_full(member_id,:)*IDF_weighting_matrix*band_adoption_full(friend_id,:)';
    temp_deno1 = sqrt(band_adoption_full(member_id,:)*IDF_weighting_matrix*band_adoption_full(member_id,:)');
    temp_deno2 = sqrt(band_adoption_full(friend_id,:)*IDF_weighting_matrix*band_adoption_full(friend_id,:)');
    sim_mf(i) = temp_nume/(temp_deno1*temp_deno2);
    sim_mf_random(i) = cosine_similarity_by_random_TF_IDF(band_adoption_full(member_id,:),band_adoption_full(friend_id,:),IDF_weighting_matrix);
end
sim_mf_rm_random = sim_mf./sim_mf_random;

sim_mn = zeros(size(neighbourlist,1),1);
for i = 1:size(neighbourlist,1)
    member_id = neighbourlist(i,1);
    neighbour_id = neighbourlist(i,2);
    temp_nume = band_adoption_full(member_id,:)*IDF_weighting_matrix*band_adoption_full(neighbour_id,:)';
    temp_deno1 = sqrt(band_adoption_full(member_id,:)*IDF_weighting_matrix*band_adoption_full(member_id,:)');
    temp_deno2 = sqrt(band_adoption_full(neighbour_id,:)*IDF_weighting_matrix*band_adoption_full(neighbour_id,:)');
    sim_mn(i) = temp_nume/(temp_deno1*temp_deno2);
end

%%
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

% MN_MF_combined = zeros(length(members_in_neighbourlist),2);
avg_sim_mn = [];
avg_sim_mf = [];
avg_sim_mf_rm_random = [];
sum(ismember(members_in_neighbourlist,members_in_friendlist))     % all members_in_neighbourlist appear in members_in_friendlist
for j=1:length(members_in_neighbourlist)
    member_index = find(members_in_neighbourlist(j)==members_in_friendlist);
    member_neighbour_rows = members_row_start_neighbourlist(j):members_row_end_neighbourlist(j);
    member_friend_rows = members_row_start_friendlist(member_index):members_row_end_friendlist(member_index);
    avg_sim_mn = [avg_sim_mn mean(sim_mn(member_neighbour_rows))];
    avg_sim_mf = [avg_sim_mf mean(sim_mf(member_friend_rows))];
    avg_sim_mf_rm_random = [avg_sim_mf_rm_random mean(sim_mf_rm_random(member_friend_rows))];
%     avg_shared_MN = mean(X_MN(member_neighbour_rows).*BP_MN(member_neighbour_rows));
%     avg_shared_MF = mean(X_MN(member_neighbour_rows).*BP_MN(member_neighbour_rows));
end

% avg_sim_mf = avg_sim_mf';

%%
avg_sim_nf = [];
avg_sim_nf_rm_random = [];
% NF_member_index = zeros(500000,1);
% sim_nf_collect = zeros(500000,1);
% count_ind = 0;
for n=1:length(members_in_neighbourlist)
    member_index = find(members_in_neighbourlist(n)==members_in_friendlist);
    member_neighbour_rows = members_row_start_neighbourlist(n):members_row_end_neighbourlist(n);
    member_friend_rows = members_row_start_friendlist(member_index):members_row_end_friendlist(member_index);
    neighbours_ind = neighbourlist(member_neighbour_rows,2);
    friends_ind = friendlist(member_friend_rows,2);
    per_member_similarity_nf = [];
    per_member_similarity_nf_rm_random = [];
    for i =1:length(neighbours_ind)
        neighbour_vec = band_adoption_full(neighbours_ind(i),:);
        for j = 1:length(friends_ind)
            friend_vec = band_adoption_full(friends_ind(j),:);
            sim_nf_scalar = cosine_similarity_TF_IDF(neighbour_vec,friend_vec,IDF_weighting_matrix);
            sim_nf_random_scalar = cosine_similarity_by_random_TF_IDF(neighbour_vec,friend_vec,IDF_weighting_matrix);
            sim_nf_random_rm_scalar = sim_nf_scalar/sim_nf_random_scalar;
            per_member_similarity_nf = [per_member_similarity_nf sim_nf_scalar];
            per_member_similarity_nf_rm_random = [per_member_similarity_nf_rm_random sim_nf_random_rm_scalar];
%             count_ind = count_ind+1;
%             sim_nf_collect(count_ind) = sim_nf_scalar;
%             NF_member_index(count_ind) = member_index(1);
        end
    end
    avg_sim_nf = [avg_sim_nf mean(per_member_similarity_nf)];
    avg_sim_nf_rm_random = [avg_sim_nf_rm_random mean(per_member_similarity_nf_rm_random)];
end

% NF_member_index = NF_member_index(1:count_ind);
% sim_nf_collect = sim_nf_collect(1:count_ind);

compareSim_mf_nf = avg_sim_mf'-avg_sim_nf';
sum(compareSim_mf_nf)

mean(avg_sim_mf)
mean(avg_sim_nf)

mean(avg_sim_mf_rm_random)
mean(avg_sim_nf_rm_random)

%% testing
% run two paired sample t test
[h,p] = ttest(avg_sim_mf,avg_sim_nf)
[h,p] = ttest(avg_sim_mf,avg_sim_nf,'Alpha',0.01)

[h,p] = ttest(avg_sim_mf_rm_random,avg_sim_nf_rm_random)
[h,p] = ttest(avg_sim_mf_rm_random,avg_sim_nf_rm_random,'Alpha',0.01)

% t test assumptions
% https://statistics.laerd.com/statistical-guides/independent-t-test-statistical-guide.php

% view normality
hist(avg_sim_mf)
hist(avg_sim_nf)
hist(avg_sim_mf-avg_sim_nf)

hist(avg_sim_mf_rm_random)
hist(avg_sim_nf_rm_random)
hist(avg_sim_mf_rm_random-avg_sim_nf_rm_random)

% view standard deviations
std(avg_sim_nf)
std(avg_sim_mf)

std(avg_sim_nf_rm_random)
std(avg_sim_mf_rm_random)

% Test for Equal Variances Using Levene’s Test
avg_sim = [avg_sim_mf'; avg_sim_nf'];
avg_sim = full(avg_sim);
group = zeros(size(avg_sim));
group(1:length(avg_sim_mf)) = 1;
p = vartestn(avg_sim,group,'TestType','LeveneAbsolute')

% run one sample test (avg_sim_mf-avg_sim_nf)
[h,p] = ttest(avg_sim_mf-avg_sim_nf,0,'Alpha',0.01)

[h,p] = ttest(avg_sim_mf_rm_random-avg_sim_nf_rm_random,0,'Alpha',0.01)

% run two independent samples t-test (samples do not need to be the same length)
[h,p,ci,stats] = ttest2(avg_sim_mf,avg_sim_nf)

[h,p,ci,stats] = ttest2(avg_sim_mf_rm_random,avg_sim_nf_rm_random)

% run Wilcoxon rank sum test (equivalent to Mann-Whitney U-test)
[p,h,stats] = ranksum(avg_sim_mf,avg_sim_nf)

[p,h,stats] = ranksum(avg_sim_mf_rm_random,avg_sim_nf_rm_random) 

% % run two independent samples t-test without averaging neighbours and friends for per member (members with more friends and/or neighbours over-sampled)
% [h,p,ci,stats] = ttest2(avg_sim_mf,sim_nf_collect)
% std(sim_nf_collect)
% std(avg_sim_mf)
% mean(sim_nf_collect)
% mean(avg_sim_mf)

% w: note: similarity computed based on only 6046 bands (with more than 1 year active and 100 listens) for all individuals
% (members, neighbours and friends) not 440000+ bands from the raw dataset

