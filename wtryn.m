% 
% temp_vec_17617085 = (1:size(prep_reverse_indices,1))';
% sparse_mat_17617085 = sparse(prep_reverse_indices(:,1),prep_reverse_indices(:,2),temp_vec_17617085,8320,6046);
% max(max(sparse_mat_17617085))
% 
% mat_sparse_8320 = sparse(17617085,8320);
% none0s_matp = mat_export(mat_export(:,7)~=0,:);
% mat_sparse_8320_row = zeros(size(none0s_matp,1),1);
% mat_sparse_8320_col = zeros(size(none0s_matp,1),1);
% mat_sparse_8320_val = zeros(size(none0s_matp,1),1);
% add_bug =0;
% for r = 1:size(none0s_matp,1)
% % for r = 5457403:size(none0s_matp,1)
%     pri_row = sparse_mat_17617085(none0s_matp(r,1),none0s_matp(r,3));
%     if pri_row == 257206
%         add_bug = add_bug+1;
%         continue
%     end
%     if rep_sparse_mat(none0s_matp(r,1),none0s_matp(r,3))==0
%         mat_sparse_8320_row(r) = prep_reverse_indices(pri_row,5)+(none0s_matp(r,4)-prep_reverse_indices(pri_row,3));
%         mat_sparse_8320_col(r) = none0s_matp(r,2);
%         mat_sparse_8320_val(r) = none0s_matp(r,7);
%     else
%         ind = rep_sparse_mat(none0s_matp(r,1),none0s_matp(r,3));
%         rows = prep_reverse_indices(rep_rows_cell{ind},:);
%         row_ind = find(rows(:,4)>=none0s_matp(r,4) & rows(:,3)<=none0s_matp(r,4));
%         if ~isempty(row_ind)
%             mat_sparse_8320_row(r) = rows(row_ind,5)+(none0s_matp(r,4)-rows(row_ind,3));
%             mat_sparse_8320_col(r) = none0s_matp(r,2);
%             mat_sparse_8320_val(r) = none0s_matp(r,7);
%         end
%     end
% end
%    


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

mean(Im_17617085)
mean(Ik_17617085)
