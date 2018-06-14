% see Onenote Similarity measure literature_User-Based Neighborhood Models_C.C.Aggarwal
% use listen times as ratings, cosine based method
% combi = combntns(1:3,2);
% combi

clc
clear all

% user_dyad = combntns(1:I,2);
% save('user_dyad_8320.mat','user_dyad', '-v7.3');           % #user 8320
% load('user_dyad_8320.mat')
friendlist = csvread('new_friendlist_8088.csv',1,0);
neighbourlist = csvread('neighbourlist_6585.csv',1,0);
memberlist = unique(neighbourlist(:,1));
all_userlist = unique([neighbourlist(:); friendlist(:)]);

i_friend_neighbour_cell = cell(length(memberlist),1);
i_non_friend_non_neighbour_cell = cell(length(memberlist),1);
for i = 1:length(memberlist)
    i_friend = friendlist(find(friendlist(:,1)==memberlist(i)),2);
    i_neighbour = friendlist(find(neighbourlist(:,1)==memberlist(i)),2);
    i_friend_neighbour_cell{i} = [i_friend; i_neighbour];
    i_non_friend_non_neighbour_cell{i} = all_userlist(find(~ismember(all_userlist,[i_friend; i_neighbour; memberlist(i)])));
end

non_f_n_ind = zeros(length(memberlist),1);
for j = 1:length(memberlist)
    non_f_n_ind(j) = length(i_non_friend_non_neighbour_cell{j});
end
cum_non_f_n_ind = cumsum(non_f_n_ind);

non_friend_neighbour_list = zeros(cum_non_f_n_ind(end),2);
non_friend_neighbour_list(1:non_f_n_ind(1),2) = i_non_friend_non_neighbour_cell{1};
non_friend_neighbour_list(1:non_f_n_ind(1),1) = memberlist(1);
for k = 2:length(memberlist)
    non_friend_neighbour_list((cum_non_f_n_ind(k-1)+1):cum_non_f_n_ind(k),2) = i_non_friend_non_neighbour_cell{k};
    non_friend_neighbour_list((cum_non_f_n_ind(k-1)+1):cum_non_f_n_ind(k),1) = memberlist(k);
end

user_dyad = non_friend_neighbour_list;

%%
friend_listens = csvread('listen_times.csv');
neigh_listens = csvread('listen_neighbours.csv',1,0);
listens = [friend_listens; neigh_listens];
% listens(:,3) = 1; %

% I = 8320
I = max(listens(:,1))
J = 6046
user_band_listen_mat = sparse(listens(:,1),listens(:,2),listens(:,3),I,J);

document_frequency = sum(user_band_listen_mat>0,1)/size(user_band_listen_mat,1);
IDF = 1-log(document_frequency);
% IDF_weight_matrix = diag(IDF.^2);
IDF_weight_matrix = sparse(1:J,1:J,IDF,J,J);
% IDF_weight_matrix = 1;

tic

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

save('Cosine_similarity_scores_ind_non_f_n_listens.mat','Cosine_similarity_scores_ind', '-v7.3');
save('Cosine_similarity_scores_non_f_n_listens.mat','Cosine_similarity_scores', '-v7.3');
save('Cosine_user_uv_ind_non_f_n_listens.mat','Cosine_user_uv_ind', '-v7.3');

% load('Jaccard_user_uv_ind_8320.mat')
% load('Jaccard_similarity_scores_ind_8320.mat')
% load('Jaccard_similarity_scores_8320.mat')

hist(Cosine_similarity_scores)
mean(Cosine_similarity_scores)

% lastfm_neighbours_match_scores = csvread('neighbourlist_6585_match.csv',1,0);
% match_neighbours = lastfm_neighbours_match_scores(:,3);
% 
% corr(Cosine_similarity_scores(1:50),match_neighbours(1:50))

%% listens
% ans =

%    0.1617


% ans =

%    0.3665

%% listen 1
% ans =

%     0.1256


% ans =

%    -0.1613

