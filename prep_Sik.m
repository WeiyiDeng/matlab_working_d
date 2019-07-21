clc
clear

% from similarity_measure_friendlist8088.m
load('cosine_similarity_scores_ind_friends8088_listen1.mat');
load('cosine_similarity_scores_friends8088_listen1.mat');
load('cosine_user_uv_ind_friends8088_listen1.mat');

cosine_similarity_sparse_mat = ...
    sparse(cosine_user_uv_ind(:,1),cosine_user_uv_ind(:,2),cosine_similarity_scores_ind,8320,8320);

% load('prep_reverse_indices.mat')

user_dyad = csvread('new_friendlist_8088.csv',1,0); 

%%
load('matp_friend_reverse.mat');
predict_trend = csvread('predict_trend_log4061_lenient.csv',1,0);

index = find(ismember(matp(:,3),predict_trend(:,1)));

matp = matp(index,:);

load('dummies40.mat');
% size(dummies40)
dummies40 = dummies40(index,:);
dummies_prep = dummies40.*repmat(matp(:,8),1,4);
clearvars dummies40
dummy_prep = dummies_prep(:,1) | dummies_prep(:,2);         % 1-20

% matp = [member_p friend_p band_p timeobs DV prob_adopt_week new_week_diff A_week_ijt];

mat_export = [matp(:,1:6) dummy_prep];

load('indx_innov.mat');
load('dummy_innov_mat.mat');

%%
% store_m = [];
% store_j = [];
% store_t_start = matp(1,4);
% store_t_end = [];
% store_r_start = 1;
% store_r_end = [];
% m = matp(1,1);
% j = matp(1,3);
% for r = 2:size(matp,1)
%     prev_m = m;
%     prev_j = j;
%     m = matp(r,1);
%     j = matp(r,3);
%     if prev_j ~= j || matp(r,2)~=matp(r-1,2)
%             store_m = [store_m prev_m];
%             store_j = [store_j prev_j];
%             store_t_start = [store_t_start matp(r,4)];
%             store_t_end = [store_t_end matp(r-1,4)];
%             store_r_start = [store_r_start r];
%             store_r_end = [store_r_end r-1];
%     end
% end
% store_m = [store_m matp(r,1)];            
% store_j = [store_j matp(r,3)];
% store_t_end = [store_t_end matp(r,4)];
% store_r_end = [store_r_end r];
% 
% store_m_index = [];
% m_index = 1;
% m_prev = store_m(1);
% for m = 1:length(store_m)
%     m_curr = store_m(m);
%     if m_curr~=m_prev
%         m_index = m_index+1;
%     end
%     store_m_index = [store_m_index m_index];
%     m_prev = m_curr;
% end

% prep_reverse_indices = [store_m' store_j' store_t_start' store_t_end' store_r_start' store_r_end' store_m_index'];
% save('prep_reverse_indices_17617085.mat','prep_reverse_indices','-v7.3');
load('prep_reverse_indices_17617085.mat')

% load('matp_friend_reverse.mat');

[C,ia,ic] = unique(prep_reverse_indices(:,1:2),'rows');
size(prep_reverse_indices,1)-length(C) 
% prep_reverse_indices = prep_reverse_indices(ia,:);
rep_index = ic(~ismember(ic,ia));   % 1077
% rep_index = unique(rep_index);
rep_rows_cell = cell(length(rep_index),1);
rep_m = zeros(length(rep_index),1);
rep_j = zeros(length(rep_index),1);
rep_ind = zeros(length(rep_index),1);
for i = 1:length(rep_index)
    m = prep_reverse_indices(rep_index(i),1);
    j = prep_reverse_indices(rep_index(i),2);
    rep_rows_cell{i} = find(prep_reverse_indices(:,1)==m & prep_reverse_indices(:,2)==j);
    rep_m(i) = m;
    rep_j(i) = j;
    rep_ind(i) = i;
end
rep_mj = [rep_m rep_j];
[C2,ia2,ic2] = unique(rep_mj,'rows');
rm_index = ic2(~ismember(ic2,ia2));
rep_m(rm_index) =[];
rep_j(rm_index) =[];
rep_ind(rm_index) =[];
rep_mj2 = [rep_m rep_j];
length(unique(rep_mj2,'rows'))
[C3,ia3,ic3] = unique(rep_mj2,'rows');
rm_index2 = ic3(~ismember(ic3,ia3));
rep_m(rm_index2) =[];
rep_j(rm_index2) =[];
rep_ind(rm_index2) =[];
rep_sparse_mat = sparse(rep_m,rep_j,rep_ind,8320,6046);
% rep_sparse_mat = sparse(rep_m,rep_j,1,8320,6046); 
% [Ir,Ic,Iv] = find(rep_sparse_mat>1);

temp_vec_17617085 = (1:size(prep_reverse_indices,1))';
sparse_mat_17617085 = sparse(prep_reverse_indices(:,1),prep_reverse_indices(:,2),temp_vec_17617085,8320,6046);
max(max(sparse_mat_17617085))

none0s_matp = mat_export(mat_export(:,7)~=0,:);
mat_sparse_8088_row = zeros(size(none0s_matp,1),1);
mat_sparse_8088_col = zeros(size(none0s_matp,1),1);
mat_sparse_8088_val = zeros(size(none0s_matp,1),1);
add_bug =0;
for r = 1:size(none0s_matp,1)
% for r = 5457403:size(none0s_matp,1)
    pri_row = sparse_mat_17617085(none0s_matp(r,1),none0s_matp(r,3));
    if pri_row == 257206
        add_bug = add_bug+1;
        continue
    end
    if rep_sparse_mat(none0s_matp(r,1),none0s_matp(r,3))==0
        mat_sparse_8088_row(r) = prep_reverse_indices(pri_row,5)+(none0s_matp(r,4)-prep_reverse_indices(pri_row,3));
        mat_sparse_8088_col(r) = cosine_similarity_sparse_mat(none0s_matp(r,1), none0s_matp(r,2));
        mat_sparse_8088_val(r) = none0s_matp(r,7);
    else
        ind = rep_sparse_mat(none0s_matp(r,1),none0s_matp(r,3));
        rows = prep_reverse_indices(rep_rows_cell{ind},:);
        row_ind = find(rows(:,4)>=none0s_matp(r,4) & rows(:,3)<=none0s_matp(r,4));
        if ~isempty(row_ind)
            mat_sparse_8088_row(r) = rows(row_ind,5)+(none0s_matp(r,4)-rows(row_ind,3));
            mat_sparse_8088_col(r) = cosine_similarity_sparse_mat(none0s_matp(r,1), none0s_matp(r,2));
            mat_sparse_8088_val(r) = none0s_matp(r,7);
        end
    end
end

sum(mat_sparse_8088_row==0)
sum(mat_sparse_8088_col==0)
rm_indices = mat_sparse_8088_row==0;
mat_sparse_8088_row(rm_indices)=[];
mat_sparse_8088_col(rm_indices)=[];
mat_sparse_8088_val(rm_indices)=[];
mat_sparse_8088 = sparse(mat_sparse_8088_row,mat_sparse_8088_col,mat_sparse_8088_val,17617085,8088);

save('mat_sparse_8088.mat','mat_sparse_8088','-v7.3');

%%
innov = csvread('EAi3_lenient.csv');

m_innov = zeros(size(cosine_user_uv_ind,1),1);
f_innov = zeros(size(cosine_user_uv_ind,1),1);
for m = 1:size(cosine_user_uv_ind,1)
    m_innov(m) = innov(cosine_user_uv_ind(m,1),2);
    f_innov(m) = innov(cosine_user_uv_ind(m,2),2);
end
    
save('m_innov.mat','m_innov','-v7.3');
save('f_innov.mat','f_innov','-v7.3');
    
%% Im and mean Ik
load('matp_friend_reverse.mat');
predict_trend = csvread('predict_trend_log4061_lenient.csv',1,0);

index = find(ismember(matp(:,3),predict_trend(:,1)));
matp = matp(index,:);

load('indx.mat');
load('dummy_mat.mat');

matp = matp(indx,:);

Im_17617085 = zeros(size(matp,1),1);
for r = 1:size(matp,1)
    Im_17617085(r) = innov(matp(r,1),2);
end

friends8088dyads_sparse_mat = sparse(cosine_user_uv_ind(:,1),cosine_user_uv_ind(:,2),1,8320,8320);

Ik_17617085 = zeros(size(matp,1),1);
for r = 1:size(matp,1)
    row_m = friends8088dyads_sparse_mat(matp(r,1),:);
    k_indices = find(row_m);
    innov_ks = innov(k_indices,2);
    Ik_17617085(r) = mean(innov_ks);
end

sum(isnan(Ik_17617085))

save('Im_17617085.mat','Im_17617085','-v7.3');
save('Ik_17617085.mat','Ik_17617085','-v7.3');


