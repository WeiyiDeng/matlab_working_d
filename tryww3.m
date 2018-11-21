load('matrix_IadoptJ_week.mat');
load('number_rows_for_DV.mat');
load('similarity_m_friend_cells.mat');
% load('m_friend_adoption_multiply_similarity_cells_Tij.mat');
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

% save('X_m_friend_adoption_times_similarity_combine.mat','X_m_friend_adoption_times_similarity_combine','-v7.3');
save('X_m_friend_adoption_times_similarity_combineD2.mat','X_m_friend_adoption_times_similarity_combineD2','-v7.3');
save('X_m_friend_adoption_times_similarity_combineD3.mat','X_m_friend_adoption_times_similarity_combineD3','-v7.3');
% load('X_m_friend_adoption_times_similarity_combine.mat');
% load('X_m_friend_adoption_times_similarity_combineD2.mat');
% load('X_m_friend_adoption_times_similarity_combineD3.mat');