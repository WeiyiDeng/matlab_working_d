similarity_m_neighbour_cells = cell(size(member_start_end,1),1);
m_neighbour_adoption_multiply_similarity_cells = cell(size(member_start_end,1),1);
m_neighbour_adoption_multiply_similarity_cells_Tij = cell(size(member_start_end,1),1);
for m = 1:size(member_start_end,1)
    m_neighbours_ID = Cosine_user_uv_ind_nm(find(Cosine_user_uv_ind_nm(:,1)==member_start_end(m,1)),2);
    similarity_m_neighbour_cells{m} = Cosine_similarity_score_nm(find(Cosine_user_uv_ind_nm(:,1)==member_start_end(m,1)));
%     m_neighbour_cells{m} = cell(size(artist_start,1),1);
    m_neighbour_adoption_multiply_similarity_cells{m} = cell(size(artist_start,1),1);
    m_neighbour_adoption_multiply_similarity_cells_Tij{m} = cell(size(artist_start,1),1);
    for j = 1:size(artist_start,1)
%         m_neighbour_cells{m}{j} = nm_adopt_6046cells{j}(m_neighbours_ID,:);
        m_neighbour_cells_mj = nm_adopt_6046cells{j}(m_neighbours_ID,:);
        temp = similarity_m_neighbour_cells{m}'*m_neighbour_cells_mj;
        m_neighbour_adoption_multiply_similarity_cells{m}{j} = sparse(1,find(temp),temp(find(temp)),1,MAX_T);
        if matrix_IJ_week_start(m,j)>matrix_IJ_week_end(m,j)              % for the bands appear after member becomes inactive, do not include
            m_neighbour_adoption_multiply_similarity_cells_Tij{m}{j} = [];
        elseif matrix_IadoptJ_week(m,j)>0                                 % if not truncated
            m_neighbour_adoption_multiply_similarity_cells_Tij{m}{j} = ...
                m_neighbour_adoption_multiply_similarity_cells{m}{j}(matrix_IJ_week_start(m,j):min(matrix_IJ_week_end(m,j),matrix_IadoptJ_week(m,j)));
        else
            m_neighbour_adoption_multiply_similarity_cells_Tij{m}{j} = ...
                m_neighbour_adoption_multiply_similarity_cells{m}{j}(matrix_IJ_week_start(m,j):matrix_IJ_week_end(m,j));
        end
    end
end