band_pop = csvread('POP_counts.csv',1,0);
band_pop(:,4) = band_pop(:,4)-104;
band_pop_mat = sparse(band_pop(:,3),band_pop(:,4),band_pop(:,2),max(band_pop(:,3)),max(band_pop(:,4)));

band_pop_sum_year = sparse(max(band_pop(:,3)),max(band_pop(:,4)));
for t = 1:51
    band_pop_sum_year(:,t) = sum(band_pop_mat(:,1:t),2);
end
for t = 52:size(band_pop_sum_year,2)
    band_pop_sum_year(:,t) = sum(band_pop_mat(:,(t-51):t),2);
end
pop = zeros(size(matp,1),1);
for r = 1:length(pop)
    pop(r) = band_pop_mat_year(matp(r,3),matp(r,4));
end
