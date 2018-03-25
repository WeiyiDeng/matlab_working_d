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

friendlist = csvread('new_friendlist_8088.csv',1,0);
friendlist_temp = friendlist(end,:);
friendlist(end,:) = [];
friendlist = [friendlist_temp; friendlist];

neighbourlist = csvread('neighbourlist_6585.csv',1,0);        % 6585 member & neighbour pairs

%%
row_start1 = 1;
row_end1 = [];
temp_member_id = friendlist(1,1);
for i = 1:size(friendlist,1)
    if friendlist(i,1) == temp_member_id
    else
        row_end1 = [row_end1 i-1];
        row_start1 = [row_start1 i];
        temp_member_id = friendlist(i,1);
    end
end
row_end1 = [row_end1 i];
collect_members_id = friendlist(row_start1,1);
        
row_start2 = 1;
row_end2 = [];
temp_member_id = neighbourlist(1,1);
for i = 1:size(neighbourlist,1)
    if neighbourlist(i,1) == temp_member_id
    else
        row_end2 = [row_end2 i-1];
        row_start2 = [row_start2 i];
        temp_member_id = neighbourlist(i,1);
    end
end
row_end2 = [row_end2 i];
collect_members_id2 = neighbourlist(row_start2,1);

collect_friends_per_member = cell(size(collect_members_id));
num_friends_per_member = zeros(size(collect_members_id));
for i = 1:length(collect_members_id)
    collect_friends_per_member{i} = friendlist(row_start1(i):row_end1(i),2);
    num_friends_per_member(i) = length(friendlist(row_start1(i):row_end1(i),2));
end

collect_neighbours_per_member = cell(size(collect_members_id2));
num_neighbours_per_member = zeros(size(collect_members_id2));
for i = 1:length(collect_members_id2)
    collect_neighbours_per_member{i} = neighbourlist(row_start2(i):row_end2(i),2);
    num_neighbours_per_member(i) = length(neighbourlist(row_start2(i):row_end2(i),2));
end

% 
row_id = zeros(500000,1);
col_id = zeros(500000,1);
A_neighbours_full_T = zeros(500000,1);
ind = 0;
for i = 1:length(collect_members_id)
    [fixed_row_id,col_id,s] = find(band_adoption_full(collect_members_id(i),:));
    for j = 1:length(col_id)
        A_neighbours_full_T_temp = sum(band_adoption_full(collect_neighbours_per_member{i,1},:));
        ind = ind+1;
        A_neighbours_full_T(ind) = A_neighbours_full_T_temp;
        row_id(ind) = fixed_row_id;
        col_id(ind) = col_id;
    end
end
