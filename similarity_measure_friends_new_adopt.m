% see Onenote Similarity measure literature_User-Based Neighborhood Models_C.C.Aggarwal
% use listen times as ratings, cosine based method
% combi = combntns(1:3,2);
% combi

clc
clear all

% friend_listens = csvread('listen_times.csv');
friend_listens = csvread('test_fm_artist_id_adopt.csv',1,0);
neigh_listens = [];
listens = [friend_listens; neigh_listens];
% listens(:,3) = 1;

% I = 8320
I = max(listens(:,1))
J = max(listens(:,2))
user_band_listen_mat = sparse(listens(:,1),listens(:,2),listens(:,3),I,J);

% user_dyad = combntns(1:I,2);
% save('user_dyad_8320.mat','user_dyad', '-v7.3');           % #user 8320
% load('user_dyad_8320.mat')
user_dyad_match = csvread('test_dt2.csv',1,0);
users_in_listens = unique(listens(:,1));
temp = ismember(user_dyad_match(:,1),users_in_listens) & ismember(user_dyad_match(:,2),users_in_listens);
sum(temp)
user_dyad_match = user_dyad_match(find(temp),:);             % remove users not in listens data

user_dyad = user_dyad_match(:,1:2);

document_frequency = sum(user_band_listen_mat>0,1)/size(user_band_listen_mat,1);
IDF = 1-log(document_frequency);
% IDF_weight_matrix = diag(IDF.^2);
IDF_weight_matrix = sparse(1:J,1:J,IDF,J,J);
% IDF_weight_matrix = 1;

tic

% Pearson_correlation = zeros(size(user_dyad,1),1);
New_adoption = zeros(size(user_dyad,1),1);
% Pearson_correlation = zeros(600000,1);
% for k = 1:6000
for k = 1:size(user_dyad,1)
    u = user_dyad(k,1);                                    % index user u
    v = user_dyad(k,2);                                    % index user v
    listen_u = user_band_listen_mat(u,:);                  
    listen_v = user_band_listen_mat(v,:);
%     New_adoption(k) = sum(listen_u.*listen_v>0 & listen_u<listen_v);
    New_adoption(k) = sum(listen_u.*listen_v>0);
%     norm_listen_u = listen_u/(sum(listen_u));              
%     norm_listen_v = listen_v/(sum(listen_v));              
%     co_rated_item_ind = find(user_band_listen_mat(u,:).*user_band_listen_mat(v,:));        % two vectors !!
%     if isempty(co_rated_item_ind)==0
%     if sum(listen_u)~=0 && sum(listen_v)~=0
% %         J_similarity(k) = Jaccard_similarity(listen_u,listen_v);
%         Cosine_similarity(k) = cosine_similarity_TF_IDF(listen_u,listen_v,IDF_weight_matrix);
% %         Cosine_similarity(k) = cosine_similarity_TF_IDF(listen_u,listen_v,IDF_weight_matrix);
% %         Cosine_similarity(k) = cosine_similarity_TF_IDF(norm_listen_u,norm_listen_v,IDF_weight_matrix);      % whether or not to normalize these two vectors gives the same cosine similarity result
% %         u_rated_items = user_band_listen_mat(u,co_rated_item_ind);
% %         v_rated_items = user_band_listen_mat(v,co_rated_item_ind);
% %         u_rates_mean_centering = u_rated_items-mean(u_rated_items);
% %         v_rates_mean_centering = v_rated_items-mean(v_rated_items);
% %         Pearson_correlation(k) = u_rates_mean_centering*v_rates_mean_centering'...
% %             /(sqrt(u_rates_mean_centering*u_rates_mean_centering')*sqrt(v_rates_mean_centering*v_rates_mean_centering'));
% %         Cosine_similarity(k) = u_rated_items*v_rated_items'...
% %             /(sqrt(u_rated_items*u_rated_items')*sqrt(v_rated_items*v_rated_items'));
%     else
% %         Pearson_correlation(k) = NaN;
% %         Cosine_similarity(k) = NaN;
%         Cosine_similarity(k) = 0;
%     end
end
% similarity_scores_ind = find(~isnan(Pearson_correlation));
% similarity_scores = Pearson_correlation(similarity_scores_ind);


New_adoption_friend_scores_ind = find(~isnan(New_adoption));
New_adoption_friend_scores = New_adoption(New_adoption_friend_scores_ind);
New_adoption_friend_user_uv_ind = user_dyad(New_adoption_friend_scores_ind,:);

toc

% save('New_adoption_friend_scores_ind_all_bands.mat','New_adoption_friend_scores_ind', '-v7.3');
% save('New_adoption_friend_scores_all_bands.mat','New_adoption_friend_scores', '-v7.3');
% save('New_adoption_friend_user_uv_ind_all_bands.mat','New_adoption_friend_user_uv_ind', '-v7.3');

% load('Jaccard_user_uv_ind_8320.mat')
% load('Jaccard_similarity_scores_ind_8320.mat')
% load('Jaccard_similarity_scores_8320.mat')

hist(New_adoption_friend_scores)
mean(New_adoption_friend_scores)

% lastfm_neighbours_match_scores = csvread('neighbourlist_6585_match.csv',1,0);
% match_neighbours = lastfm_neighbours_match_scores(:,3);
% match_neighbours = user_dyad_match(:,3);
% hist(match_neighbours)

% corr(Cosine_similarity_scores, match_neighbours)

%% listen 1
% ans =

%     0.1094

%% listen
% ans =

%    0.4631