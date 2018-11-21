% see Onenote Similarity measure literature_User-Based Neighborhood Models_C.C.Aggarwal
% use listen times as ratings, cosine based method
% combi = combntns(1:3,2);
% combi

clc
clear all

% friend_listens = csvread('listen_times.csv');
friend_listens = csvread('all_test_fm_id.csv',1,0);
neigh_listens = csvread('all_test_nm_id.csv',1,0);
listens = [friend_listens; neigh_listens];
% listens(:,3) = 1;

% I = 8320
I = max(listens(:,1))
J = max(listens(:,2))
user_band_listen_mat = sparse(listens(:,1),listens(:,2),listens(:,3),I,J);

% user_dyad = combntns(1:I,2);
% save('user_dyad_8320.mat','user_dyad', '-v7.3');           % #user 8320
% load('user_dyad_8320.mat')
MFID = csvread('ALL_USERS_dt2.csv',1,0);
MNID = csvread('ALL_USERS_dt3.csv',1,0);

temp = unique(MFID(:,1));
temp2 = unique(MNID(:,1));
member_list = temp(ismember(temp,temp2));

max_user = max([MNID(:,1);MNID(:,2);MFID(:,1);MFID(:,2)]);
prep = repmat(member_list',max_user,1);
prep = prep(:);
prep2 = repmat((1:max_user)',length(member_list),1);
combi = [prep prep2];
rm_index = find(prep==prep2);
combi(rm_index,:) = [];
rm_index2 = ismember(combi,MFID,'rows');
combi(rm_index2,:) = [];
user_dyad_full = combi;

b = 300000;
n = size(user_dyad_full,1);
cells = mat2cell(user_dyad_full,diff([0:b:n-1,n]));

% users_in_listens = unique(listens(:,1));
% temp = ismember(user_dyad_match(:,1),users_in_listens) & ismember(user_dyad_match(:,2),users_in_listens);
% sum(temp)
% user_dyad_match = user_dyad_match(find(temp),:);             % remove users not in listens data
% 
% user_dyad = user_dyad_match(:,1:2);

document_frequency = sum(user_band_listen_mat>0,1)/size(user_band_listen_mat,1);
IDF = 1-log(document_frequency);
% IDF_weight_matrix = diag(IDF.^2);
IDF_weight_matrix = sparse(1:J,1:J,IDF,J,J);
% IDF_weight_matrix = 1;

tic

%%
user_dyad = cells{4};

% Pearson_correlation = zeros(size(user_dyad,1),1);
Cosine_similarity = zeros(size(user_dyad,1),1);
% Pearson_correlation = zeros(600000,1);
% for k = 1:6000
for k = 1:size(user_dyad,1)
    u = user_dyad(k,1);                                    % index user u
    v = user_dyad(k,2);                                    % index user v
    listen_u = user_band_listen_mat(u,:);                  
    listen_v = user_band_listen_mat(v,:);
%     norm_listen_u = listen_u/(sum(listen_u));              
%     norm_listen_v = listen_v/(sum(listen_v));              
%     co_rated_item_ind = find(user_band_listen_mat(u,:).*user_band_listen_mat(v,:));        % two vectors !!
%     if isempty(co_rated_item_ind)==0
    if sum(listen_u)~=0 && sum(listen_v)~=0
%         J_similarity(k) = Jaccard_similarity(listen_u,listen_v);
        Cosine_similarity(k) = cosine_similarity_TF_IDF(listen_u,listen_v,IDF_weight_matrix);
%         Cosine_similarity(k) = cosine_similarity_TF_IDF(listen_u,listen_v,IDF_weight_matrix);
%         Cosine_similarity(k) = cosine_similarity_TF_IDF(norm_listen_u,norm_listen_v,IDF_weight_matrix);      % whether or not to normalize these two vectors gives the same cosine similarity result
%         u_rated_items = user_band_listen_mat(u,co_rated_item_ind);
%         v_rated_items = user_band_listen_mat(v,co_rated_item_ind);
%         u_rates_mean_centering = u_rated_items-mean(u_rated_items);
%         v_rates_mean_centering = v_rated_items-mean(v_rated_items);
%         Pearson_correlation(k) = u_rates_mean_centering*v_rates_mean_centering'...
%             /(sqrt(u_rates_mean_centering*u_rates_mean_centering')*sqrt(v_rates_mean_centering*v_rates_mean_centering'));
%         Cosine_similarity(k) = u_rated_items*v_rated_items'...
%             /(sqrt(u_rated_items*u_rated_items')*sqrt(v_rated_items*v_rated_items'));
    else
%         Pearson_correlation(k) = NaN;
%         Cosine_similarity(k) = NaN;
        Cosine_similarity(k) = 0;
    end
end
% similarity_scores_ind = find(~isnan(Pearson_correlation));
% similarity_scores = Pearson_correlation(similarity_scores_ind);


Cosine_similarity_scores_ind = find(~isnan(Cosine_similarity));
Cosine_similarity_scores = Cosine_similarity(Cosine_similarity_scores_ind);
Cosine_user_uv_ind = user_dyad(Cosine_similarity_scores_ind,:);

toc

%%
save('Cosine_similarity_ind_job4.mat','Cosine_similarity_scores_ind', '-v7.3');
save('Cosine_similarity_scores_job4.mat','Cosine_similarity_scores', '-v7.3');
save('Cosine_user_uv_ind_job4.mat','Cosine_user_uv_ind', '-v7.3');

% load('Cosine_similarity_scores_ind_neighbours_new_listens_TFIDF.mat')
% load('Cosine_similarity_scores_neighbours_new_listens_TFIDF.mat')
% load('Cosine_user_uv_ind_neighbours_new_listens_TFIDF.mat')

% hist(Cosine_similarity_scores)
% mean(Cosine_similarity_scores)
% 
% % lastfm_neighbours_match_scores = csvread('neighbourlist_6585_match.csv',1,0);
% % match_neighbours = lastfm_neighbours_match_scores(:,3);
% match_neighbours = user_dyad_match(:,3);
% hist(match_neighbours)
% 
% corr(Cosine_similarity_scores, match_neighbours)

%% listen
% ans =

%    0.7570


% ans =

%  -2.9259e-04

%% listen 1
% ans =

%    0.0563


% ans =

%    0.2044