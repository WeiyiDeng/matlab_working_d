clc
clear

load('count_shared_non_nfbands.mat')
% load('Cosine_similarity_scores_non_f_n_listens.mat')
load('Cosine_similarity_scores_non_f_n_listen1.mat')

index_non_nf = find(~isnan(Cosine_similarity_scores) & ~isnan(count_shared_non_nfbands));
beta_non_nf = Cosine_similarity_scores(index_non_nf)\count_shared_non_nfbands(index_non_nf)