load('ALL_USERS_NM.mat')
load('ALL_USERS_NM_similarity_scores.mat')
load('m_r_start.mat')
load('m_r_end.mat')
load('N_id_cells.mat')
load('memberlistNM.mat')
ALL_NEW_OLD_ID = csvread('ALL_USER_NEW_OLD8320_ID.csv',1,0);
load('prep_reverse_indices.mat')

similarities_prep_for_matp = [];
for m = 1:size(N_id_cells,1)
    temp = prep_reverse_indices(find(prep_reverse_indices(:,7)==m),1);
    m_old = temp(1);
    m_new = ALL_NEW_OLD_ID(m_old,2);
    m_N_list = N_id_cells{m};
    if ~isempty(m_N_list)
        m_N_similarities = zeros(length(m_N_list),1);
        for n = 1:length(m_N_list)
            m_N_index = find(ALL_USERS_NM(:,1)==m_new & ALL_USERS_NM(:,2)==m_N_list(n));
            m_N_similarities(n) = ALL_USERS_NM_similarity_scores(m_N_index);
        end
    else
        m_N_similarities = [];
    end
    similarities_prep_for_matp = [similarities_prep_for_matp; m_N_similarities];
end

save('similarities_prep_for_matp.mat','similarities_prep_for_matp','-v7.3');    
