clc
clear

adoptions = csvread('bandadoptions_lenient_adopt.csv',1,0);
friendlist = csvread('new_friendlist_8088.csv',1,0);

adoptions_neighbour = csvread('bandadoptions_neighbour.csv',1,0);
neighbourlist = csvread('neighbourlist_6585.csv',1,0);

adoptions_full = [adoptions(:,1:3);adoptions_neighbour];      % add neighbour adoptions to below

I = max(adoptions_full(:,1))
J = max(adoptions_full(:,2))
user_band_adoption_mat = sparse(adoptions_full(:,1),adoptions_full(:,2),adoptions_full(:,3),I,J);

count_shared_fbands = zeros(size(friendlist,1),1);
temp = 0;
temp2 = 0;
for k = 1:size(friendlist,1)
    u = friendlist(k,1);
    v = friendlist(k,2);
    u_adopt = user_band_adoption_mat(u,:);
    v_adopt = user_band_adoption_mat(v,:);
    shared_bands = u_adopt.*v_adopt > 0; % & v_adopt-u_adopt >= 0 & v_adopt-u_adopt <= 52;
    count_shared_fbands(k) = sum(shared_bands);
    temp = temp+sum(u_adopt>0);
    temp2 = temp2+sum(v_adopt>0);
end
temp/size(friendlist,1)
temp2/size(friendlist,1)

count_shared_nbands = zeros(size(neighbourlist,1),1);
temp3 = zeros(size(neighbourlist,1),1);
temp4 = zeros(size(neighbourlist,1),1);
for n = 1:size(neighbourlist,1)
    u = neighbourlist(n,1);
    v = neighbourlist(n,2);
    u_adopt = user_band_adoption_mat(u,:);
    v_adopt = user_band_adoption_mat(v,:);
    shared_bands = u_adopt.*v_adopt > 0; % & v_adopt-u_adopt >= 0 & v_adopt-u_adopt <= 52;
    count_shared_nbands(n) = sum(shared_bands);
    temp3(n) = sum(u_adopt>0);
    temp4(n) = sum(v_adopt>0);
end
mean(temp3)
mean(temp4)
   
mean(count_shared_fbands)
mean(count_shared_nbands)

friend_translist = csvread('new_old_flist.csv',1,0);
neighbour_translist = csvread('new_old_nlist.csv',1,0);

% load('Cosine_similarity_scores_neighbour_new_artist_listen1_TFIDF.mat')
% load('Cosine_user_uv_ind_neighbour_new_artist_listen1_TFIDF.mat');
% new_neighbourlist = Cosine_user_uv_ind;
% cosine_similarity_score_nm = Cosine_similarity_scores;
% 
% load('Cosine_similarity_scores_friend_new_artist_listen1_TFIDF.mat')
% load('Cosine_user_uv_ind_friend_new_artist_listen1_TFIDF.mat');
% new_friendlist = Cosine_user_uv_ind;
% cosine_similarity_score_fm = Cosine_similarity_scores;

load('Cosine_similarity_scores_new_neighbours_withTFIDF.mat')
load('Cosine_user_uv_ind_new_neighbours_withTFIDF.mat');
new_neighbourlist = Cosine_user_uv_ind;
cosine_similarity_score_nm = Cosine_similarity_scores;

load('Cosine_similarity_scores_new_friends_withTFIDF.mat')
load('Cosine_user_uv_ind_new_friends_withTFIDF.mat');
new_friendlist = Cosine_user_uv_ind;
cosine_similarity_score_fm = Cosine_similarity_scores;

mystore_n = zeros(max(neighbour_translist(:,2)),1);
mystore_n(neighbour_translist(:,2)) = neighbour_translist(:,1);
new_neighbourlist_transformed = zeros(size(new_neighbourlist));
for i = 1:size(new_neighbourlist,1)
    new_neighbourlist_transformed(i,1) = mystore_n(new_neighbourlist(i,1));
    new_neighbourlist_transformed(i,2) = mystore_n(new_neighbourlist(i,2));
end

mystore_f = zeros(max(friend_translist(:,2)),1);
mystore_f(friend_translist(:,2)) = friend_translist(:,1);
new_friendlist_transformed = zeros(size(new_friendlist));
for i = 1:size(new_friendlist,1)
    new_friendlist_transformed(i,1) = mystore_f(new_friendlist(i,1));
    new_friendlist_transformed(i,2) = mystore_f(new_friendlist(i,2));
end

similarity_nm = zeros(size(neighbourlist,1),1);
for i = 1:size(neighbourlist,1)
    search_index = find(new_neighbourlist_transformed(:,1)==neighbourlist(i,1) & new_neighbourlist_transformed(:,2)==neighbourlist(i,2));
    if isempty(search_index)
        similarity_nm(i) = NaN;
    else
        similarity_nm(i) = cosine_similarity_score_nm(search_index(1));
    end
end
 
similarity_fm = zeros(size(friendlist,1),1);
for i = 1:size(friendlist,1)
    search_index = find(new_friendlist_transformed(:,1)==friendlist(i,1) & new_friendlist_transformed(:,2)==friendlist(i,2));
    if isempty(search_index)
        similarity_fm(i) = NaN;
    else
        similarity_fm(i) = cosine_similarity_score_fm(search_index(1));
    end
end
   
%% test
length(unique(friendlist(:,1)))
length(unique(neighbourlist(:,1)))
sth = ismember(unique(friendlist(:,1)),unique(neighbourlist(:,1)));          % there is one less member in neighbourlist than in friendlist

%% new tracks & artists
uni_nlist = unique(neighbourlist(:,1));
nlist_index = cell(size(uni_nlist));
flist_index = cell(size(uni_nlist));
for i = 1:length(uni_nlist)
    n_index = find(neighbourlist(:,1)==uni_nlist(i));
    f_index = find(friendlist(:,1)==uni_nlist(i));
    nlist_index{i} = n_index;
    flist_index{i} = f_index;
end

metric = zeros(size(uni_nlist));
for d = 1:length(uni_nlist)
    n_ind = nlist_index{d};
    f_ind = flist_index{d};
    An = count_shared_nbands(n_ind);
    Sn = similarity_nm(n_ind);
    Af = count_shared_fbands(f_ind);
    Sf = similarity_fm(f_ind);
    beta_n = Sn\An;
    An_hat = beta_n.*Sf;
    metric(d) = mean(Af./An_hat);
end
    
mean(metric(~isnan(metric)))    
    
% %% 6046 artists    
% load('Cosine_similarity_scores_neighbours6585_listens.mat')
% cosine_similarity_score_6046nm = Cosine_similarity_scores;
% 
% load('cosine_similarity_scores_friends8088_listens.mat')
% cosine_similarity_score_6046fm = cosine_similarity_scores;   
% 
% uni_nlist = unique(neighbourlist(:,1));
% nlist_row = zeros(length(uni_nlist),2);
% flist_row = zeros(length(uni_nlist),2);
% for i = 1:length(uni_nlist)
%     n_rows = find(neighbourlist(:,1)==uni_nlist(i));
%     f_rows = find(friendlist(:,1)==uni_nlist(i));
%     if isempty(f_rows) && ~isempty(n_rows)
%         flist_row(i,:) = [NaN NaN];
%         nlist_row(i,:) = [min(n_rows) max(n_rows)];
%     elseif isempty(n_rows) && ~isempty(f_rows)
%         nlist_row(i,:) = [NaN NaN];
%         flist_row(i,:) = [min(f_rows) max(f_rows)];
%     elseif isempty(n_rows) && isempty(f_rows)
%         flist_row(i,:) = [NaN NaN];
%         nlist_row(i,:) = [NaN NaN];
%     else
%         nlist_row(i,:) = [min(n_rows) max(n_rows)];
%         flist_row(i,:) = [min(f_rows) max(f_rows)];
%     end
% end
% 
% max(max(nlist_row))
% max(max(flist_row))
% 
% metric6046 = zeros(size(uni_nlist));
% sth = 0;
% for d = 1:length(uni_nlist)
%     n_ind = nlist_index{d};
%     f_ind = flist_index{d};
%     if isnan(nlist_row(d,:)+flist_row(d,:))
%         metric6046(d) = NaN;
%     elseif d == 55                 % last line of friendlist8088, MEMBER_ID = 2898 and FRIEND_ID = 6613
%         An = count_shared_nbands(n_ind);
%         Sn = cosine_similarity_score_6046nm(nlist_row(d,1):nlist_row(d,2));
%         Af = count_shared_fbands(f_ind);
%         Sf = cosine_similarity_score_6046fm([1:9 8088]);
%         beta_n = Sn\An;
%         An_hat = beta_n.*Sf;
%         metric6046(d) = nanmean(Af./An_hat);
%     else
%         An = count_shared_nbands(n_ind);
%         Sn = cosine_similarity_score_6046nm(nlist_row(d,1):nlist_row(d,2));
%         Af = count_shared_fbands(f_ind);
%         Sf = cosine_similarity_score_6046fm(flist_row(d,1):flist_row(d,2));
%         beta_n = Sn\An;
%         An_hat = beta_n.*Sf;
%         metric6046(d) = nanmean(Af./An_hat);
%     end
% end
%     
% mean(metric6046(~isnan(metric6046)))  
%     