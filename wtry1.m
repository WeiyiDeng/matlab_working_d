mt_verif = zeros(length(mtind),1);
mtind = find(matp(:,5));
for i = 1:length(mtind)
    mt_temp = band_adopt_mat(matp(mtind(i),1),matp(mtind(i),3));
    mt_verif(i) = mt_temp==matp(mtind(i),4);
    if i>4000
        warning(msg)
    end
end