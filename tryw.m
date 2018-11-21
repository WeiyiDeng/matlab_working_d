% get the adoption times similarity with friends for each member M (M cells {J cells{all M's friends I*all time period Tij}})
% m_neighbour_cells = cell(size(member_start_end,1),1);
similarity_m_friend_cells = cell(size(member_start_end,1),1);
m_friend_adoption_multiply_similarity_cells = cell(size(member_start_end,1),1);
m_friend_adoption_multiply_similarity_cellsD1 = cell(size(member_start_end,1),1);
m_friend_adoption_multiply_similarity_cellsD4 = cell(size(member_start_end,1),1);
m_friend_adoption_multiply_similarity_cells_Tij = cell(size(member_start_end,1),1);
m_friend_adoption_multiply_similarity_cells_TijD1 = cell(size(member_start_end,1),1);
m_friend_adoption_multiply_similarity_cells_TijD4 = cell(size(member_start_end,1),1);
for m = 1:size(member_start_end,1)
    m_friend_ID = Cosine_user_uv_ind_fm(find(Cosine_user_uv_ind_fm(:,1)==member_start_end(m,2)),2);
    similarity_m_friend_cells{m} = Cosine_similarity_score_fm(find(Cosine_user_uv_ind_fm(:,1)==member_start_end(m,2)));
%     m_neighbour_cells{m} = cell(size(artist_start,1),1);
    m_friend_adoption_multiply_similarity_cells{m} = cell(size(artist_start,1),1);
    m_friend_adoption_multiply_similarity_cellsD1{m} = cell(size(artist_start,1),1);
    m_friend_adoption_multiply_similarity_cellsD4{m} = cell(size(artist_start,1),1);
    m_friend_adoption_multiply_similarity_cells_Tij{m} = cell(size(artist_start,1),1);
    m_friend_adoption_multiply_similarity_cells_TijD1{m} = cell(size(artist_start,1),1);
    m_friend_adoption_multiply_similarity_cells_TijD4{m} = cell(size(artist_start,1),1);
    for j = 1:size(artist_start,1)
%         m_neighbour_cells{m}{j} = nm_adopt_6046cells{j}(m_neighbours_ID,:);
        m_friend_cells_mj = fm_adopt_6046cells{j}(m_friend_ID,:);
        temp = similarity_m_friend_cells{m}'*m_friend_cells_mj;
        m_friend_adoption_multiply_similarity_cells{m}{j} = sparse(1,find(temp),temp(find(temp)),1,MAX_T);
        m_friend_adoption_multiply_similarity_cellsD1{m}{j} = [sparse(0,[],[],1,1) m_friend_adoption_multiply_similarity_cells{m}{j}(1:end-1)];
        m_friend_adoption_multiply_similarity_cellsD4{m}{j} = [sparse(0,[],[],1,4) m_friend_adoption_multiply_similarity_cells{m}{j}(1:end-4)];
        if matrix_IJ_week_start(m,j)>matrix_IJ_week_end(m,j)              % for the bands appear after member becomes inactive, do not include
            m_friend_adoption_multiply_similarity_cells_Tij{m}{j} = [];
            m_friend_adoption_multiply_similarity_cells_TijD1{m}{j} = [];
            m_friend_adoption_multiply_similarity_cells_TijD4{m}{j} = [];
        elseif matrix_IadoptJ_week(m,j)>0                                 % if not truncated
            m_friend_adoption_multiply_similarity_cells_Tij{m}{j} = ...
                m_friend_adoption_multiply_similarity_cells{m}{j}(matrix_IJ_week_start(m,j):min(matrix_IJ_week_end(m,j),matrix_IadoptJ_week(m,j)));
            m_friend_adoption_multiply_similarity_cells_TijD1{m}{j} = ...
                m_friend_adoption_multiply_similarity_cellsD1{m}{j}(matrix_IJ_week_start(m,j):min(matrix_IJ_week_end(m,j),matrix_IadoptJ_week(m,j)));
            m_friend_adoption_multiply_similarity_cells_TijD4{m}{j} = ...
                m_friend_adoption_multiply_similarity_cellsD4{m}{j}(matrix_IJ_week_start(m,j):min(matrix_IJ_week_end(m,j),matrix_IadoptJ_week(m,j)));
        else
            m_friend_adoption_multiply_similarity_cells_Tij{m}{j} = ...
                m_friend_adoption_multiply_similarity_cells{m}{j}(matrix_IJ_week_start(m,j):matrix_IJ_week_end(m,j));
            m_friend_adoption_multiply_similarity_cells_TijD1{m}{j} = ...
                m_friend_adoption_multiply_similarity_cellsD1{m}{j}(matrix_IJ_week_start(m,j):matrix_IJ_week_end(m,j));
            m_friend_adoption_multiply_similarity_cells_TijD4{m}{j} = ...
                m_friend_adoption_multiply_similarity_cellsD4{m}{j}(matrix_IJ_week_start(m,j):matrix_IJ_week_end(m,j));
        end
    end
end

sum_ww = 0;
for i = 1:165
    ww = 0;
    for j = 1:6046
        ww = ww + length(find(m_friend_adoption_multiply_similarity_cells_TijD4{i}{j}));
%         if ww>0
%             disp(j)
%         end
    end
    sum_ww = sum_ww+ww;
end
disp(sum_ww)

