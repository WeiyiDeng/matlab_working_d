similarity_m_friend_cells = cell(size(member_start_end,1),1);
m_friend_adoption_multiply_similarity_cells = cell(size(member_start_end,1),1);
m_friend_adoption_multiply_similarity_cellsD2 = cell(size(member_start_end,1),1);
m_friend_adoption_multiply_similarity_cellsD3 = cell(size(member_start_end,1),1);
m_friend_adoption_multiply_similarity_cells_Tij = cell(size(member_start_end,1),1);
m_friend_adoption_multiply_similarity_cells_TijD2 = cell(size(member_start_end,1),1);
m_friend_adoption_multiply_similarity_cells_TijD3 = cell(size(member_start_end,1),1);
for m = 1:size(member_start_end,1)
    m_friend_ID = Cosine_user_uv_ind_fm(find(Cosine_user_uv_ind_fm(:,1)==member_start_end(m,2)),2);
    similarity_m_friend_cells{m} = Cosine_similarity_score_fm(find(Cosine_user_uv_ind_fm(:,1)==member_start_end(m,2)));
%     m_neighbour_cells{m} = cell(size(artist_start,1),1);
    m_friend_adoption_multiply_similarity_cells{m} = cell(size(artist_start,1),1);
    m_friend_adoption_multiply_similarity_cellsD2{m} = cell(size(artist_start,1),1);
    m_friend_adoption_multiply_similarity_cellsD3{m} = cell(size(artist_start,1),1);
    m_friend_adoption_multiply_similarity_cells_Tij{m} = cell(size(artist_start,1),1);
    m_friend_adoption_multiply_similarity_cells_TijD2{m} = cell(size(artist_start,1),1);
    m_friend_adoption_multiply_similarity_cells_TijD3{m} = cell(size(artist_start,1),1);
    for j = 1:size(artist_start,1)
%         m_neighbour_cells{m}{j} = nm_adopt_6046cells{j}(m_neighbours_ID,:);
        m_friend_cells_mj = fm_adopt_6046cells{j}(m_friend_ID,:);
        temp = similarity_m_friend_cells{m}'*m_friend_cells_mj;
        m_friend_adoption_multiply_similarity_cells{m}{j} = sparse(1,find(temp),temp(find(temp)),1,MAX_T);
        m_friend_adoption_multiply_similarity_cellsD2{m}{j} = [sparse(0,[],[],1,2) m_friend_adoption_multiply_similarity_cells{m}{j}(1:end-2)];         % adopted 2 week ago
        m_friend_adoption_multiply_similarity_cellsD3{m}{j} = [sparse(0,[],[],1,3) m_friend_adoption_multiply_similarity_cells{m}{j}(1:end-3)];         % adopted 3 week ago
        if matrix_IJ_week_start(m,j)>matrix_IJ_week_end(m,j)              % for the bands appear after member becomes inactive, do not include
            m_friend_adoption_multiply_similarity_cells_Tij{m}{j} = [];
            m_friend_adoption_multiply_similarity_cells_TijD2{m}{j} = [];
            m_friend_adoption_multiply_similarity_cells_TijD3{m}{j} = [];
        elseif matrix_IadoptJ_week(m,j)>0                                 % if not truncated
            m_friend_adoption_multiply_similarity_cells_Tij{m}{j} = ...
                m_friend_adoption_multiply_similarity_cells{m}{j}(matrix_IJ_week_start(m,j):min(matrix_IJ_week_end(m,j),matrix_IadoptJ_week(m,j)));
            m_friend_adoption_multiply_similarity_cells_TijD2{m}{j} = ...
                m_friend_adoption_multiply_similarity_cellsD2{m}{j}(matrix_IJ_week_start(m,j):min(matrix_IJ_week_end(m,j),matrix_IadoptJ_week(m,j)));
            m_friend_adoption_multiply_similarity_cells_TijD3{m}{j} = ...
                m_friend_adoption_multiply_similarity_cellsD3{m}{j}(matrix_IJ_week_start(m,j):min(matrix_IJ_week_end(m,j),matrix_IadoptJ_week(m,j)));
        else
            m_friend_adoption_multiply_similarity_cells_Tij{m}{j} = ...
                m_friend_adoption_multiply_similarity_cells{m}{j}(matrix_IJ_week_start(m,j):matrix_IJ_week_end(m,j));
            m_friend_adoption_multiply_similarity_cells_TijD2{m}{j} = ...
                m_friend_adoption_multiply_similarity_cellsD2{m}{j}(matrix_IJ_week_start(m,j):matrix_IJ_week_end(m,j));
            m_friend_adoption_multiply_similarity_cells_TijD3{m}{j} = ...
                m_friend_adoption_multiply_similarity_cellsD3{m}{j}(matrix_IJ_week_start(m,j):matrix_IJ_week_end(m,j));
        end
    end
end

% save('similarity_m_friend_cells.mat','similarity_m_friend_cells','-v7.3');
% save('m_friend_adoption_multiply_similarity_cells.mat','m_friend_adoption_multiply_similarity_cells','-v7.3');
% save('m_friend_adoption_multiply_similarity_cells_Tij.mat','m_friend_adoption_multiply_similarity_cells_Tij','-v7.3');
save('m_friend_adoption_multiply_similarity_cellsD2.mat','m_friend_adoption_multiply_similarity_cellsD2','-v7.3');
save('m_friend_adoption_multiply_similarity_cells_TijD2.mat','m_friend_adoption_multiply_similarity_cells_TijD2','-v7.3');
save('m_friend_adoption_multiply_similarity_cellsD3.mat','m_friend_adoption_multiply_similarity_cellsD3','-v7.3');
save('m_friend_adoption_multiply_similarity_cells_TijD3.mat','m_friend_adoption_multiply_similarity_cells_TijD3','-v7.3');

clear similarity_m_friend_cells
clear m_friend_adoption_multiply_similarity_cells
clear m_friend_adoption_multiply_similarity_cells_Tij
clear m_friend_adoption_multiply_similarity_cellsD2
clear m_friend_adoption_multiply_similarity_cells_TijD2
clear m_friend_adoption_multiply_similarity_cellsD3
clear m_friend_adoption_multiply_similarity_cells_TijD3

load('similarity_m_friend_cells.mat');
load('m_friend_adoption_multiply_similarity_cells_Tij.mat');
load('m_friend_adoption_multiply_similarity_cells_TijD2.mat');
load('m_friend_adoption_multiply_similarity_cells_TijD3.mat');

% X_m_friend_adoption_times_similarity_combine = [];
X_m_friend_adoption_times_similarity_combineD2 = [];
X_m_friend_adoption_times_similarity_combineD3 = [];
for m = 1:size(matrix_IadoptJ_week,1)
    for j = 1:size(matrix_IadoptJ_week,2)
        if number_rows_for_DV(m,j)~= 0
%             X_m_friend_adoption_times_similarity_sparse = m_friend_adoption_multiply_similarity_cells_Tij{m}{j}';
            X_m_friend_adoption_times_similarity_sparseD2 = m_friend_adoption_multiply_similarity_cells_TijD2{m}{j}';
            X_m_friend_adoption_times_similarity_sparseD3 = m_friend_adoption_multiply_similarity_cells_TijD3{m}{j}';
%             X_m_friend_adoption_times_similarity_combine = ...
%                 [X_m_friend_adoption_times_similarity_combine; X_m_friend_adoption_times_similarity_sparse];
            X_m_friend_adoption_times_similarity_combineD2 = ...
                [X_m_friend_adoption_times_similarity_combineD2; X_m_friend_adoption_times_similarity_sparseD2];
            X_m_friend_adoption_times_similarity_combineD3 = ...
                [X_m_friend_adoption_times_similarity_combineD3; X_m_friend_adoption_times_similarity_sparseD3];
        else
        end
    end
end

save('X_m_friend_adoption_times_similarity_combine.mat','X_m_friend_adoption_times_similarity_combine','-v7.3');
save('X_m_friend_adoption_times_similarity_combineD2.mat','X_m_friend_adoption_times_similarity_combineD2','-v7.3');
save('X_m_friend_adoption_times_similarity_combineD3.mat','X_m_friend_adoption_times_similarity_combineD3','-v7.3');
% load('X_m_friend_adoption_times_similarity_combine.mat');
% load('X_m_friend_adoption_times_similarity_combineD2.mat');
% load('X_m_friend_adoption_times_similarity_combineD3.mat');
