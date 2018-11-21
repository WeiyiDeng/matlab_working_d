All_Cosine_similarity_scores = [];
All_Cosine_similarity_scores_ind = [];
All_Cosine_user_uv_ind = [];
for i = 1:11
    name = sprintf('Cosine_similarity_scores_job%d.mat',i);
    name2 = sprintf('Cosine_similarity_ind_job%d.mat',i);
    name3 = sprintf('Cosine_user_uv_ind_job%d.mat',i);
    load(name)
    load(name2)
    load(name3)
    All_Cosine_similarity_scores = [All_Cosine_similarity_scores; Cosine_similarity_scores];
    All_Cosine_similarity_scores_ind = [All_Cosine_similarity_scores_ind; Cosine_similarity_scores_ind];
    All_Cosine_user_uv_ind = [All_Cosine_user_uv_ind; Cosine_user_uv_ind];
end
save('All_Cosine_similarity_scores.mat','All_Cosine_similarity_scores','-v7.3');
save('All_Cosine_similarity_ind.mat','All_Cosine_similarity_ind','-v7.3');
save('All_Cosine_user_uv_ind.mat','All_Cosine_user_uv_ind','-v7.3');

memberlist = unique(All_Cosine_user_uv_ind(:,1));
m_start = 1;
m_end = [];
prev_m = All_Cosine_user_uv_ind(1,1);
for i = 2:size(All_Cosine_user_uv_ind,1)
    m = All_Cosine_user_uv_ind(i,1);
    if m~=prev_m
        m_start = [m_start; i];
        m_end = [m_end; i-1];
    end
    prev_m = m;
end
m_end = [m_end; i];

K = 5;                                              % top K most similar neighbours
top_neighbours_row = zeros(length(memberlist),K);
top_similarities = zeros(length(memberlist),K);
top_neighbours = zeros(length(memberlist),K);
for m = 1:length(memberlist)
% for m = 1    
    temp = All_Cosine_similarity_scores(m_start(m):m_end(m));
    [temp_ind,originalpos] = sort(temp, 'descend' );
    n = temp_ind(1:K);
    p = originalpos(1:K)+m_start(m)-1;
    top_neighbours_row(m,:) = p;
    top_neighbours(m,:) = All_Cosine_user_uv_ind(p,2);
    top_similarities(m,:) = n;
end
    
ALL_USERS_NM = csvread('ALL_USERS_dt3.csv',1,0);
ALL_USERS_NM = ALL_USERS_NM(ismember(ALL_USERS_NM(:,1),memberlist),:);
memberlistNM = unique(ALL_USERS_NM(:,1));
sum(memberlistNM == memberlist)

ALL_USERS_NM_similarity_scores = zeros(length(ALL_USERS_NM),1);
for k = 1:length(ALL_USERS_NM)
    index = find(All_Cosine_user_uv_ind(:,1)==ALL_USERS_NM(k,1) &...
        All_Cosine_user_uv_ind(:,2)==ALL_USERS_NM(k,2));
    if isempty(index)
        ALL_USERS_NM_similarity_scores(k) = NaN;
    else
        ALL_USERS_NM_similarity_scores(k) = All_Cosine_similarity_scores(index);
    end
end

rm_indices = isnan(ALL_USERS_NM_similarity_scores);
ALL_USERS_NM_similarity_scores(rm_indices)=[];
ALL_USERS_NM(rm_indices,:) = []; 

m_r_start = 1;
m_r_end = [];
m_id = ALL_USERS_NM(1,1);
prev_m = ALL_USERS_NM(1,1);
for i = 2:size(ALL_USERS_NM,1)
    m = ALL_USERS_NM(i,1);
    if m~=prev_m
        m_r_start = [m_r_start; i];
        m_id = [m_id; ALL_USERS_NM(i,1)];
        m_r_end = [m_r_end; i-1];
    end
    prev_m = m;
end
m_r_end = [m_r_end; i];

memberlistNM = m_id;
% save('memberlistNM.mat','memberlistNM','-v7.3');
% save('m_r_start.mat','m_r_start','-v7.3'); 
% save('m_r_end.mat','m_r_end','-v7.3'); 
% save('ALL_USERS_NM.mat','ALL_USERS_NM','-v7.3'); 
% save('ALL_USERS_NM_similarity_scores.mat','ALL_USERS_NM_similarity_scores','-v7.3'); 

load('memberlistNM.mat')
load('m_r_start.mat')
load('m_r_end.mat')
load('ALL_USERS_NM.mat')
load('ALL_USERS_NM_similarity_scores.mat')



