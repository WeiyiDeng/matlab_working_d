clc
clear

load('matp_friend_reverse.mat');
adoptions = csvread('bandadoptions_lenient_adopt.csv',1,0);
length(unique(matp(:,1)))
length(unique(matp(:,2)))
ww = unique(matp(:,2));
friendlist = csvread('new_friendlist_8088.csv',1,0);
w = unique(matp(:,1));
sum(ismember(w,friendlist(:,1)))
length(unique(friendlist(:,1)))
sum(ismember(ww,friendlist(:,2)))
length(unique(friendlist(:,2)))

load('memberlistNM.mat')
load('ALL_USERS_NM.mat')
load('ALL_USERS_NM_similarity_scores.mat')
load('m_r_start.mat')                             
load('m_r_end.mat')

ALL_NEW_OLD_ID = csvread('ALL_USER_NEW_OLD8320_ID.csv',1,0);
ALL_USERS_NM = csvread('ALL_USERS_dt3.csv',1,0);
ALL_USERS_FM = csvread('ALL_USERS_dt2.csv',1,0);

load('matp_friend_reverse.mat');
all_adoptions = csvread('dt7_all_usersID6046_WEEKID.csv',1,0);
max(all_adoptions(:,2))
max(all_adoptions(:,1))
max(max(ALL_USERS_FM))
max(max(ALL_USERS_NM))
[C, ia, ic] = unique(all_adoptions(:,1:2),'rows');
all_adoptions_no_duplicate = all_adoptions(ia,:);

sum(ismember(unique(ALL_USERS_NM(:,1)),unique(ALL_USERS_FM(:,1))))

u1 = unique(ALL_USERS_NM(:,1));
u2 = unique(ALL_USERS_FM(:,1));
sth=u1(ismember(u1,u2));
sum(ismember(u1,u2))

new_id_members = ALL_NEW_OLD_ID(ismember(ALL_NEW_OLD_ID(:,1),unique(matp(:,1))),2);

ALL_USERS_NM_143 = ALL_USERS_NM(ismember(ALL_USERS_NM(:,1),new_id_members),:);
length(unique(ALL_USERS_NM_143(:,1)))
ALL_USERS_FM_143 = ALL_USERS_FM(ismember(ALL_USERS_FM(:,1),new_id_members),:);
length(unique(ALL_USERS_FM_143(:,1)))

all_adoptions_sparse_mat = sparse(all_adoptions_no_duplicate(:,1), all_adoptions_no_duplicate(:,2), all_adoptions_no_duplicate(:,3),...
    max(all_adoptions_no_duplicate(:,1)),max(all_adoptions_no_duplicate(:,2)));
max(max(all_adoptions_sparse_mat))

count_shared_week_differences_NM = cell(size(ALL_USERS_NM_143,1),1);
for n = 1:size(ALL_USERS_NM_143,1)
    member_id = ALL_USERS_NM_143(n,1);
    neighbour_id = ALL_USERS_NM_143(n,2);
    m_adopt_vec = all_adoptions_sparse_mat(member_id,:);
    n_adopt_vec = all_adoptions_sparse_mat(neighbour_id,:);
    shared_adopt_index = find(m_adopt_vec.*n_adopt_vec);
    m_adopt_shared = m_adopt_vec(shared_adopt_index);
    n_adopt_shared = n_adopt_vec(shared_adopt_index);
    count_shared_week_differences_NM{n} = m_adopt_shared-n_adopt_shared;
end
count_shared_week_differences_vec_NM = [];
for i = 1:size(count_shared_week_differences_NM)
    temp = count_shared_week_differences_NM{i};
    within_40_weeks = temp(abs(temp)<=40);
    count_shared_week_differences_vec_NM = [count_shared_week_differences_vec_NM within_40_weeks];
end
figure
hist(count_shared_week_differences_vec_NM)

count_shared_week_differences_FM = cell(size(ALL_USERS_FM_143,1),1);
for f = 1:size(ALL_USERS_FM_143,1)
    member_id = ALL_USERS_FM_143(f,1);
    friend_id = ALL_USERS_FM_143(f,2);
    m_adopt_vec = all_adoptions_sparse_mat(member_id,:);
    f_adopt_vec = all_adoptions_sparse_mat(friend_id,:);
    shared_adopt_index = find(m_adopt_vec.*f_adopt_vec);
    m_adopt_shared = m_adopt_vec(shared_adopt_index);
    f_adopt_shared = f_adopt_vec(shared_adopt_index);
    count_shared_week_differences_FM{f} = m_adopt_shared-f_adopt_shared;
end
count_shared_week_differences_vec_FM = [];
for i = 1:size(count_shared_week_differences_FM)
    temp = count_shared_week_differences_FM{i};
    within_40_weeks = temp(abs(temp)<=40);
    count_shared_week_differences_vec_FM = [count_shared_week_differences_vec_FM within_40_weeks];
end
figure
hist(count_shared_week_differences_vec_FM)

csvwrite('count_shared_week_differences_vec_FM.csv', full(count_shared_week_differences_vec_FM))   % plot 1
csvwrite('count_shared_week_differences_vec_NM.csv', full(count_shared_week_differences_vec_NM))   % plot 3

load('All_Cosine_user_uv_ind.mat')
load('All_Cosine_similarity_scores.mat')
load('ALL_USERS_NM_similarity_scores.mat')

sum(ismember(ALL_USERS_FM(:,1),All_Cosine_user_uv_ind(:)))
sum(ismember(ALL_USERS_NM(:,1),All_Cosine_user_uv_ind(:)))

all_similarities_sparse_mat = sparse(All_Cosine_user_uv_ind(:,1), All_Cosine_user_uv_ind(:,2), ...
    All_Cosine_similarity_scores, max(All_Cosine_user_uv_ind(:,1)), max(All_Cosine_user_uv_ind(:,2)));

similarities_143_NM = zeros(size(ALL_USERS_NM_143,1),1);
for n = 1:size(ALL_USERS_NM_143,1)
    similarities_143_NM(n) = all_similarities_sparse_mat(ALL_USERS_NM_143(n,1),ALL_USERS_NM_143(n,2));
end

similarities_143_FM = zeros(size(ALL_USERS_FM_143,1),1);
for f = 1:size(ALL_USERS_FM_143,1)
    similarities_143_FM(f) = all_similarities_sparse_mat(ALL_USERS_FM_143(f,1),ALL_USERS_FM_143(f,2));
end

count_shared_adoptions_per_member_within40weeks_NM = zeros(size(ALL_USERS_NM_143,1),1);
for i = 1:size(count_shared_week_differences_NM)
    temp = count_shared_week_differences_NM{i};
    count_shared_adoptions_per_member_within40weeks_NM(i) = sum(abs(temp)<=40);
end
count_shared_adoptions_per_member_within40weeks_FM = zeros(size(ALL_USERS_FM_143,1),1);
for i = 1:size(count_shared_week_differences_FM)
    temp = count_shared_week_differences_FM{i};
    count_shared_adoptions_per_member_within40weeks_FM(i) = sum(abs(temp)<=40);
end

count_shared_adoptions_per_member_within20weeks_NM = zeros(size(ALL_USERS_NM_143,1),1);
for i = 1:size(count_shared_week_differences_NM)
    temp = count_shared_week_differences_NM{i};
    count_shared_adoptions_per_member_within20weeks_NM(i) = sum(abs(temp)<=20);
end
count_shared_adoptions_per_member_within20weeks_FM = zeros(size(ALL_USERS_FM_143,1),1);
for i = 1:size(count_shared_week_differences_FM)
    temp = count_shared_week_differences_FM{i};
    count_shared_adoptions_per_member_within20weeks_FM(i) = sum(abs(temp)<=20);
end

count_shared_adoptions_per_member_within12weeks_NM = zeros(size(ALL_USERS_NM_143,1),1);
for i = 1:size(count_shared_week_differences_NM)
    temp = count_shared_week_differences_NM{i};
    count_shared_adoptions_per_member_within12weeks_NM(i) = sum(abs(temp)<=12);
end
count_shared_adoptions_per_member_within12weeks_FM = zeros(size(ALL_USERS_FM_143,1),1);
for i = 1:size(count_shared_week_differences_FM)
    temp = count_shared_week_differences_FM{i};
    count_shared_adoptions_per_member_within12weeks_FM(i) = sum(abs(temp)<=12);
end

count_shared_adoptions_per_member_within4weeks_NM = zeros(size(ALL_USERS_NM_143,1),1);
for i = 1:size(count_shared_week_differences_NM)
    temp = count_shared_week_differences_NM{i};
    count_shared_adoptions_per_member_within4weeks_NM(i) = sum(abs(temp)<=4);
end
count_shared_adoptions_per_member_within4weeks_FM = zeros(size(ALL_USERS_FM_143,1),1);
for i = 1:size(count_shared_week_differences_FM)
    temp = count_shared_week_differences_FM{i};
    count_shared_adoptions_per_member_within4weeks_FM(i) = sum(abs(temp)<=4);
end

scatter(similarities_143_NM, count_shared_adoptions_per_member_within40weeks_NM);
scatter(similarities_143_FM, count_shared_adoptions_per_member_within40weeks_FM);

flist_to_ALL_USERS_ID = csvread('test_flist_to_ALL_USERS_ID.csv',1,0);

load('Cosine_similarity_scores_friend_new_listens_TFIDF.mat')
load('Cosine_user_uv_ind_friend_new_listens_TFIDF.mat')

similarities_143_FM_real = zeros(size(ALL_USERS_FM_143,1),1);
for m = 1: size(ALL_USERS_FM_143,1)
    member_id = flist_to_ALL_USERS_ID(find(flist_to_ALL_USERS_ID(:,2)==ALL_USERS_FM_143(m,1)),1);
    friend_id = flist_to_ALL_USERS_ID(find(flist_to_ALL_USERS_ID(:,2)==ALL_USERS_FM_143(m,2)),1);
    if ~isempty(friend_id) && ~isempty(member_id)
        ind = find(Cosine_user_uv_ind(:,1)==member_id & Cosine_user_uv_ind(:,2)==friend_id);
        similarities_143_FM_real(m) = Cosine_similarity_scores(ind(1));
    else
        similarities_143_FM_real(m) = 0;
    end
end
scatter(similarities_143_FM_real, count_shared_adoptions_per_member_within40weeks_FM);

% plot 2
csvwrite('similarities_143_NM.csv', similarities_143_NM)
csvwrite('count_shared_adoptions_per_member_within40weeks_NM.csv', count_shared_adoptions_per_member_within40weeks_NM)
csvwrite('similarities_143_FM_real.csv', similarities_143_FM_real)
csvwrite('count_shared_adoptions_per_member_within40weeks_FM.csv', count_shared_adoptions_per_member_within40weeks_FM)


scatter([similarities_143_FM_real; similarities_143_NM],...
    [count_shared_adoptions_per_member_within40weeks_FM; count_shared_adoptions_per_member_within40weeks_NM]);

similarities_143_NM_sparse_mat = sparse(ALL_USERS_NM_143(:,1), ALL_USERS_NM_143(:,2), similarities_143_NM);
count_shared_adoptions_per_member_within20weeks_NM_sparse_mat = ...
    sparse(ALL_USERS_NM_143(:,1), ALL_USERS_NM_143(:,2), count_shared_adoptions_per_member_within20weeks_NM);

NM_member_list = unique(ALL_USERS_NM_143(:,1));
n1_shared_adoptions = zeros(size(NM_member_list));
n2_shared_adoptions = zeros(size(NM_member_list));
n1_similarity = zeros(size(NM_member_list));
n2_similarity = zeros(size(NM_member_list));
for i = 1:length(NM_member_list)
    m_id = NM_member_list(i);
    n_sim_vec = similarities_143_NM_sparse_mat(m_id,:);
    [max1, ind1] = max(n_sim_vec);
    n_sim_vec(ind1) = -Inf;
    [max2, ind2] = max(n_sim_vec);
    n1_shared_adoptions(i) = count_shared_adoptions_per_member_within20weeks_NM_sparse_mat(m_id, ind1);
    n2_shared_adoptions(i) = count_shared_adoptions_per_member_within20weeks_NM_sparse_mat(m_id, ind2);
    n1_similarity(i) = max1;
    n2_similarity(i) = max2;    
end

scatter(n1_similarity, n1_shared_adoptions);
scatter(n2_similarity, n2_shared_adoptions);  

temp = n1_shared_adoptions+n2_shared_adoptions;
plot(n1_similarity, n2_similarity, '.')
    
% https://www.mathworks.com/matlabcentral/fileexchange/5105-making-surface-plots-from-scatter-data

tri = delaunay(n1_similarity, n2_similarity);
h = trisurf(tri, n1_similarity, n2_similarity, temp);
axis vis3d

% another arrangement
n1_sim = [n1_similarity; zeros(size(n2_similarity))];
n2_sim = [zeros(size(n1_similarity)); n2_similarity];
n1_n2_shared_adopt = [n1_shared_adoptions; n2_shared_adoptions];
plot(n1_sim, n2_sim, '.')
tri = delaunay(n1_sim, n2_sim);
h = trisurf(tri, n1_sim, n2_sim, n1_n2_shared_adopt,'FaceAlpha',0.5);
axis vis3d
h.EdgeColor = 'none';

c = linspace(1,10,length(n1_sim));
scatter3(n1_sim,n2_sim,n1_n2_shared_adopt,10,c,'filled');
for i = 1:length(n1_similarity)
    col1 = [n1_sim(i);n1_sim(length(n1_similarity)+i)];
    col2 = [n2_sim(i);n2_sim(length(n1_similarity)+i)];
    col3 = [n1_n2_shared_adopt(i);n1_n2_shared_adopt(length(n1_similarity)+i)];
    line(col1,col2,col3)
end


csvwrite('n1_similarity.csv', n1_similarity)
csvwrite('n2_similarity.csv', n2_similarity)
csvwrite('n1_shared_adoptions_within20weeks.csv', n1_shared_adoptions)
csvwrite('n2_shared_adoptions_within20weeks.csv', n2_shared_adoptions)
