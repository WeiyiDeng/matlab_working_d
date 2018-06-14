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

%%
length(unique(friendlist(:,1)))
length(unique(neighbourlist(:,1)))

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

count_shared_non_nfbands = zeros(size(non_friend_neighbour_list,1),1);
% temp = 0;
% temp2 = 0;
for n = 1:size(non_friend_neighbour_list,1)
    u = non_friend_neighbour_list(n,1);
    v = non_friend_neighbour_list(n,2);
    u_adopt = user_band_adoption_mat(u,:);
    v_adopt = user_band_adoption_mat(v,:);
    shared_bands = u_adopt.*v_adopt > 0 & v_adopt-u_adopt >= 0; % & v_adopt-u_adopt <= 52;
    count_shared_non_nfbands(n) = sum(shared_bands);
%     temp = temp+sum(u_adopt>0);
%     temp2 = temp2+sum(v_adopt>0);
end
% temp/size(friendlist,1)
% temp2/size(friendlist,1)

mean(count_shared_non_nfbands)

% save('count_shared_non_nfbands.mat','count_shared_non_nfbands', '-v7.3');
load('count_shared_non_nfbands.mat')


