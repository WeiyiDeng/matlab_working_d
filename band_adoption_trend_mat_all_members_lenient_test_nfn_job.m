clc
clear all

% adoptions = csvread('bandadoptions3.csv');               % note that adoptions has not subtracted 104
adoptions = csvread('bandadoptions_lenient_adopt.csv',1,0);    
adoptions = adoptions(:,1:3);
adopt_ones = ones(size(adoptions, 1),1);

% total number of weeks: length(105:527) = 423
T = 423
I = 8320                                                 % 165 members
J = 6046
old_bandt = 104
band_adoption = cell(T,1);

for t = 1:T
    [r c v] = find(adoptions(:,3)==(t+old_bandt));    
    adopt_ind = adoptions(r,:);
    band_adoption{t} = sparse(adopt_ind(:,1),adopt_ind(:,2),adopt_ones(r),I,J);
end

band_adopt_mat = sparse(adoptions(:,1),adoptions(:,2),adoptions(:,3)-old_bandt,I,J);              

%%
all_adoptions = csvread('dt7_all_usersID6046_WEEKID.csv',1,0);

all_adoptions(:,3) = all_adoptions(:,3)-104;
J_bands = max(all_adoptions(:,2));
I_users = max(all_adoptions(:,1));
T_weeks = max(all_adoptions(:,3));

% all_adopt_cells = cell(1,J_bands);
% for j = 1:J_bands
%     j_rows = all_adoptions(find(all_adoptions(:,2)==j),:);
%     all_adopt_cells{j} = sparse(j_rows(:,1),j_rows(:,3),1,I_users,T_weeks);
% end
% [r,c,v] = find(all_adopt_cells{1});
% for j = 1:J_bands
%     all_adopt_cells{j} = spones(all_adopt_cells{j});               % set all replicated adoptions to one
% end

% save('all_adopt_cells.mat','all_adopt_cells','-v7.3');
load('all_adopt_cells.mat')

ALL_NEW_OLD_ID = csvread('ALL_USER_NEW_OLD8320_ID.csv',1,0);

friendlist = csvread('new_friendlist_8088.csv',1,0);
member_old_id = unique(friendlist(:,1));
member_new_id = ALL_NEW_OLD_ID(member_old_id,2);

% ALL_USERS_dt3
ALL_USERS_NM = csvread('ALL_USERS_dt3.csv',1,0);
ALL_USERS_FM = csvread('ALL_USERS_dt2.csv',1,0);

% load('matp_friend_reverse.mat')
% store_m = [];
% store_j = [];
% store_t_start = matp(1,4);
% store_t_end = [];
% store_r_start = 1;
% store_r_end = [];
% m = matp(1,1);
% j = matp(1,3);
% for r = 2:size(matp,1)
%     prev_m = m;
%     prev_j = j;
%     m = matp(r,1);
%     j = matp(r,3);
%     if prev_j ~= j || matp(r,2)~=matp(r-1,2)
%             store_m = [store_m prev_m];
%             store_j = [store_j prev_j];
%             store_t_start = [store_t_start matp(r,4)];
%             store_t_end = [store_t_end matp(r-1,4)];
%             store_r_start = [store_r_start r];
%             store_r_end = [store_r_end r-1];
%     end
% end
% store_m = [store_m matp(r,1)];            
% store_j = [store_j matp(r,3)];
% store_t_end = [store_t_end matp(r,4)];
% store_r_end = [store_r_end r];
% 
% store_m_index = [];
% m_index = 1;
% m_prev = store_m(1);
% for m = 1:length(store_m)
%     m_curr = store_m(m);
%     if m_curr~=m_prev
%         m_index = m_index+1;
%     end
%     store_m_index = [store_m_index m_index];
%     m_prev = m_curr;
% end

% prep_reverse_indices = [store_m' store_j' store_t_start' store_t_end' store_r_start' store_r_end' store_m_index'];
% save('prep_reverse_indices.mat','prep_reverse_indices','-v7.3');
load('prep_reverse_indices.mat')

% num_members = length(unique(prep_reverse_indices(:,7)));
% 
% member_cube_cells = cell(num_members,1);
% for m = 1:num_members
%     member_all_adopt_cells = cell(J_bands,1);
%     for j = 1:J_bands
%         member_all_adopt_cells{j} = sparse(I_users,T_weeks);
%     end
%     member_cube_cells{m} = member_all_adopt_cells;
% end
% save('member_cube_cells.mat','member_cube_cells','-v7.3');

% load('member_cube_cells.mat')
load('memberlistNM.mat')
load('ALL_USERS_NM.mat')
load('ALL_USERS_NM_similarity_scores.mat')
load('m_r_start.mat')                              % 177 m ???
load('m_r_end.mat')

% test_store = [];
% not_exist_m = [];
% N_id_cells = cell(max(prep_reverse_indices(i,7)),1);
% N_similarity_cells = cell(max(prep_reverse_indices(i,7)),1);
% for i = 1:size(prep_reverse_indices,1)
%     m = prep_reverse_indices(i,7);
%     m_new = ALL_NEW_OLD_ID(prep_reverse_indices(i,1),2);
%     if ismember(m_new,memberlistNM)
%         m_new_row = find(memberlistNM==m_new);
%         N_id_list = ALL_USERS_NM(m_r_start(m_new_row):m_r_end(m_new_row),2);
%         N_id_cells{m} = N_id_list;
%         N_similarity_cells{m} = ALL_USERS_NM_similarity_scores(m_r_start(m_new_row):m_r_end(m_new_row));
%         j = prep_reverse_indices(i,2);
%         week_interval = prep_reverse_indices(i,3):prep_reverse_indices(i,4);
%         temp = all_adopt_cells{j}(:,week_interval);
%         member_cube_cells{m}{j} = temp(N_id_list,:);
%         [row, col] = find(member_cube_cells{m}{j});
%         if ~isempty(row)
%             test_store = [test_store; length(unique(row))];
%         end
%         col = flip(col);
%         row = flip(row);
%         for r = 1:length(row)
%             member_cube_cells{m}{j}(row(r),:) = abs(week_interval-col(r)-week_interval(1));
%         end
%     else
%         not_exist_m = [not_exist_m m];
%     end
% end       
% not_exist_m = unique(not_exist_m);
% save('member_cube_cells.mat','member_cube_cells','-v7.3');
% save('N_id_cells.mat','N_id_cells','-v7.3');
% save('N_similarity_cells.mat','N_similarity_cells','-v7.3');
% save('not_exist_m.mat','not_exist_m','-v7.3');

% count_N = [];
% for i = 1:size(N_id_cells,1)
%     count_N = [count_N length(N_id_cells{i})];
% end
% max(count_N)
% sum(count_N)
% 
% count_N_end = cumsum(count_N);
% count_N_start = [1 cumsum(count_N)+1];
% count_N_start(end) = [];

% save('count_N.mat','count_N','-v7.3');
% save('count_N_start.mat','count_N_start','-v7.3');
% save('count_N_end.mat','count_N_end','-v7.3');

load('count_N.mat')
load('count_N_start.mat')
load('count_N_end.mat')
load('N_id_cells.mat')
load('N_similarity_cells.mat')
load('not_exist_m.mat')

count_N = [];
for i = 1:size(N_id_cells,1)
    count_N = [count_N length(N_id_cells{i})];
end
max(count_N)
sum(count_N)

count_N_end = cumsum(count_N);
count_N_start = [1 cumsum(count_N)+1];
count_N_start(end) = [];

% N_prep_for_matp = sparse(prep_reverse_indices(end,6),sum(count_N));
% for i = 1:size(prep_reverse_indices,1)
%     m = prep_reverse_indices(i,7);
%     j = prep_reverse_indices(i,2);
%     row_interval = prep_reverse_indices(i,5):prep_reverse_indices(i,6);
%     N_prep_for_matp(row_interval,count_N_start(m):count_N_end(m))...
%         = (member_cube_cells{m}{j})';
% end
    
clearvars member_cube_cells N_prep_for_matp
test_store = [];
not_exist_m = [];
N_id_cells = cell(max(prep_reverse_indices(:,7)),1);
N_similarity_cells = cell(max(prep_reverse_indices(:,7)),1);
N_prep_for_matp = sparse(prep_reverse_indices(end,6),6222);           % sum(count_N)
for i = 1:size(prep_reverse_indices,1)
    m = prep_reverse_indices(i,7);
    m_new = ALL_NEW_OLD_ID(prep_reverse_indices(i,1),2);
    if ismember(m_new,memberlistNM)
        m_new_row = find(memberlistNM==m_new);
        N_id_list = ALL_USERS_NM(m_r_start(m_new_row):m_r_end(m_new_row),2);
        N_id_cells{m} = N_id_list;
        N_similarity_cells{m} = ALL_USERS_NM_similarity_scores(m_r_start(m_new_row):m_r_end(m_new_row));
        j = prep_reverse_indices(i,2);
        week_interval = prep_reverse_indices(i,3):prep_reverse_indices(i,4);
        temp = all_adopt_cells{j}(:,week_interval);
        member_cube_cells{m}{j} = temp(N_id_list,:);
        [row, col] = find(member_cube_cells{m}{j});
        if ~isempty(row)
            test_store = [test_store; length(unique(row))];
        end
        col = flip(col);
        row = flip(row);
        for r = 1:length(row)
            member_cube_cells{m}{j}(row(r),:) = abs(week_interval-col(r)-week_interval(1));
        end
        row_interval = prep_reverse_indices(i,5):prep_reverse_indices(i,6);
        N_prep_for_matp(row_interval,count_N_start(m):count_N_end(m))...
            = (member_cube_cells{m}{j})';
    else
        not_exist_m = [not_exist_m m];
    end
end       
not_exist_m = unique(not_exist_m);

save('N_prep_for_matp.mat','N_prep_for_matp','-v7.3');
save('member_cube_cells.mat','member_cube_cells','-v7.3');
% save('N_id_cells.mat','N_id_cells','-v7.3');
% save('N_similarity_cells.mat','N_similarity_cells','-v7.3');
% save('not_exist_m.mat','not_exist_m','-v7.3');
