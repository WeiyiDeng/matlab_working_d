APOP_ready_vec = zeros(size(m_t_agg,1),1);
for r = 1:size(m_t_agg,1)
    APOP_ready_vec(r) = member_agg_bands_APOP_cells{member_nummer(m_t_agg(r,1))}(m_t_agg(r,2));
end
    
member_id_rec = m_t_agg(1,1);
m_end_rec = [];
m_start_rec = 1;
m_list_rec = [];
for i = 2:size(m_t_agg,1)
    prev_m = member_id_rec;
    member_id_rec = m_t_agg(i,1);
    if prev_m ~=member_id_rec
        m_end_rec = [m_end_rec; i-1];
        m_start_rec = [m_start_rec; i];
        m_list_rec = [m_list_rec; prev_m];
    else
    end
end
m_list_rec = [m_list_rec; member_id_rec];    
m_end_rec = [m_end_rec; i];

min(m_end_rec-m_start_rec)

for t = 1:10
    APOP_ready_vec_rec = []
    lag_y_rec = []
    for i = 1:length(m_list_rec)
        lag_y_temp = y_pred(m_start_rec(i):m_end_rec(i)-t);
        APOP_temp = APOP_ready_vec(m_start_rec(i)+t:m_end_rec(i),:);
        APOP_ready_vec_rec = [APOP_ready_vec_rec; APOP_temp];
        lag_y_rec = [lag_y_rec; lag_y_temp];
    end
    SimuX_rec = [ones(size(lag_y_rec,1),1) zeros(size(lag_y_rec,1),1) APOP_ready_vec_rec lag_y_rec zeros(size(lag_y_rec,1),1) APOP_ready_vec_rec.*lag_y_rec];
    y_rec = SimuX_rec*modelresults;
end

