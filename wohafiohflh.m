
diffusion_year = zeros(T,J);
for i = 1:7
    ind = (i-1)*52+1;
    year_adopt = sum(diffusion_jt(ind:ind+51,:),1);
    diffusion_year(ind:ind+51,:) = repmat(year_adopt,52,1);
end
year_adopt_8th = sum(diffusion_jt(365:423,:),1);
diffusion_year(365:423,:) = repmat(year_adopt_8th,59,1);