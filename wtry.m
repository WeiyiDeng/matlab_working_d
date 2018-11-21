clearvars member_cube_cells N_prep_for_matp
test_store = [];
not_exist_m = [];
N_id_cells = cell(max(prep_reverse_indices(:,7)),1);
N_similarity_cells = cell(max(prep_reverse_indices(:,7)),1);
N_prep_for_matp = sparse(prep_reverse_indices(end,6),6222);           % sum(count_N)
for i = 1:size(prep_reverse_indices,1)
% for i = 165517:165518
% for i = 1:1
    m = prep_reverse_indices(i,7);
    m_new = ALL_NEW_OLD_ID(prep_reverse_indices(i,1),2);
    if ismember(m_new,memberlistNM)
        m_new_row = find(memberlistNM==m_new);
        N_id_list = ALL_USERS_NM(m_r_start(m_new_row):m_r_end(m_new_row),2);
        N_id_cells{m} = N_id_list;
        N_similarity_cells{m} = ALL_USERS_NM_similarity_scores(m_r_start(m_new_row):m_r_end(m_new_row));
        j = prep_reverse_indices(i,2);
        week_interval = prep_reverse_indices(i,3):prep_reverse_indices(i,4);
        temp = all_adopt_cells{j}(:,week_interval);
        member_cube_cells{m}{j} = temp(N_id_list,:);
        [row, col] = find(member_cube_cells{m}{j});
        if ~isempty(row)
            test_store = [test_store; length(unique(row))];
        end
        col = flip(col);
        row = flip(row);
        for r = 1:length(row)
            member_cube_cells{m}{j}(row(r),:) = abs(week_interval-col(r)-week_interval(1));
        end
        row_interval = prep_reverse_indices(i,5):prep_reverse_indices(i,6);
        N_prep_for_matp(row_interval,count_N_start(m):count_N_end(m))...
            = (member_cube_cells{m}{j})';
    else
        not_exist_m = [not_exist_m m];
    end
end       
not_exist_m = unique(not_exist_m);