paper3_DV = zeros(size(m_t_agg,1),1);
for i = 1:size(m_t_agg,1)
    ind = find(m_t_usage(:,1)==m_t_agg(i,1) & m_t_usage(:,2)==m_t_agg(i,2));
    if length(ind)>1
        disp(i)
    elseif isempty(ind)
        paper3_DV(i) = 0;
    else
        paper3_DV(i) = m_t_usage(ind,3);
    end
end