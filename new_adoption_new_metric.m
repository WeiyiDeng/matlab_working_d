% see Onenote Similarity measure literature_User-Based Neighborhood Models_C.C.Aggarwal
% use listen times as ratings, cosine based method
% combi = combntns(1:3,2);
% combi

clc
clear all

%%
% load('New_adoption_neighbour_scores_all_bands.mat')
% load('New_adoption_neighbour_user_uv_ind_all_bands.mat');
% 
% load('New_adoption_friend_scores_all_bands.mat')
% load('New_adoption_friend_user_uv_ind_all_bands.mat');

% load('New_adoption_neighbour_scores_all_bands_52.mat')
% load('New_adoption_neighbour_user_uv_ind_all_bands_52.mat');
% 
% load('New_adoption_friend_scores_all_bands_52.mat')
% load('New_adoption_friend_user_uv_ind_all_bands_52.mat');
 
% load('New_adoption_neighbour_scores_all_bands_12.mat')
% load('New_adoption_neighbour_user_uv_ind_all_bands_4.mat');
% 
% load('New_adoption_friend_scores_all_bands_12.mat')
% load('New_adoption_friend_user_uv_ind_all_bands_4.mat');

load('New_adoption_neighbour_scores_all_bands_4.mat')
load('New_adoption_neighbour_user_uv_ind_all_bands_4.mat');

load('New_adoption_friend_scores_all_bands_4.mat')
load('New_adoption_friend_user_uv_ind_all_bands_4.mat');

%%

% load('Cosine_similarity_scores_neighbour_new_artist_listens_TFIDF.mat')
% load('Cosine_user_uv_ind_neighbour_new_artist_listens_TFIDF.mat');
% new_neighbourlist = Cosine_user_uv_ind;
% cosine_similarity_score_nm = Cosine_similarity_scores;
% 
% load('Cosine_similarity_scores_friend_new_artist_listens_TFIDF.mat')
% load('Cosine_user_uv_ind_friend_new_artist_listens_TFIDF.mat');
% new_friendlist = Cosine_user_uv_ind;
% cosine_similarity_score_fm = Cosine_similarity_scores;

% load('Cosine_similarity_scores_neighbour_new_artist_listen1_TFIDF.mat')
% load('Cosine_user_uv_ind_neighbour_new_artist_listen1_TFIDF.mat');
% new_neighbourlist = Cosine_user_uv_ind;
% cosine_similarity_score_nm = Cosine_similarity_scores;
% 
% load('Cosine_similarity_scores_friend_new_artist_listen1_TFIDF.mat')
% load('Cosine_user_uv_ind_friend_new_artist_listen1_TFIDF.mat');
% new_friendlist = Cosine_user_uv_ind;
% cosine_similarity_score_fm = Cosine_similarity_scores;

% load('Cosine_similarity_scores_neighbours_new_listens_TFIDF.mat')           % nanmean(Sn_vec) < nanmean(Sf_vec)
% load('Cosine_user_uv_ind_neighbours_new_listens_TFIDF.mat');
% new_neighbourlist = Cosine_user_uv_ind;
% cosine_similarity_score_nm = Cosine_similarity_scores;
% 
% load('Cosine_similarity_scores_friend_new_listens_TFIDF.mat')
% load('Cosine_user_uv_ind_friend_new_listens_TFIDF.mat');
% new_friendlist = Cosine_user_uv_ind;
% cosine_similarity_score_fm = Cosine_similarity_scores;

load('Cosine_similarity_scores_neighbour_new_listen1_TFIDF.mat')
load('Cosine_user_uv_ind_neighbours_new_listens_TFIDF.mat');
new_neighbourlist = Cosine_user_uv_ind;
cosine_similarity_score_nm = Cosine_similarity_scores;

load('Cosine_similarity_scores_friend_new_listen1_TFIDF.mat')
load('Cosine_user_uv_ind_friend_new_listen1_TFIDF.mat');
new_friendlist = Cosine_user_uv_ind;
cosine_similarity_score_fm = Cosine_similarity_scores;

index_f = find(~isnan(cosine_similarity_score_fm) & ~isnan(New_adoption_friend_scores));
index_n = find(~isnan(cosine_similarity_score_nm) & ~isnan(New_adoption_neighbour_scores));

disp(['corr_f :    ' num2str(corr(cosine_similarity_score_fm(index_f),New_adoption_friend_scores(index_f))) ''])
disp(['corr_n :    ' num2str(corr(cosine_similarity_score_nm(index_n),New_adoption_neighbour_scores(index_n))) ''])



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