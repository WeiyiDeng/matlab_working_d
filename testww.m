clear

load('Cosine_similarity_scores_neighbour_new_artist_listen1_TFIDF.mat')
cosine_similarity_score_nm_artist_lis1 = Cosine_similarity_scores;

load('Cosine_similarity_scores_friend_new_artist_listen1_TFIDF.mat')
cosine_similarity_score_fm_artist_lis1 = Cosine_similarity_scores;

load('Cosine_similarity_scores_neighbours_new_listens_TFIDF.mat')
cosine_similarity_score_nm_listens = Cosine_similarity_scores;

load('Cosine_similarity_scores_friend_new_listens_TFIDF.mat')
cosine_similarity_score_fm_listens = Cosine_similarity_scores;

load('Cosine_similarity_scores_neighbour_new_listen1_TFIDF.mat');
cosine_similarity_score_nm_listen1 = Cosine_similarity_scores;

load('Cosine_similarity_scores_friend_new_listen1_TFIDF.mat');
cosine_similarity_score_fm_listen1 = Cosine_similarity_scores;

corr(cosine_similarity_score_nm_artist_lis1,cosine_similarity_score_nm_listens)         % r = 0.2558
corr(cosine_similarity_score_nm_artist_lis1,cosine_similarity_score_nm_listen1)         % r = 0.7592
corr(cosine_similarity_score_nm_listen1,cosine_similarity_score_nm_listens)             % r = 0.2903

corr(cosine_similarity_score_fm_artist_lis1,cosine_similarity_score_fm_listens)         % r = 0.4419
corr(cosine_similarity_score_fm_artist_lis1,cosine_similarity_score_fm_listen1)         % r = 0.8429
corr(cosine_similarity_score_fm_listen1,cosine_similarity_score_fm_listens)             % r = 0.4139