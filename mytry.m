DV_IJ_sparse_combine = [];
for i = 1:size(position_IadoptJ_in_rows,1)
    for j = 1:size(position_IadoptJ_in_rows,2)
        if number_rows_for_DV(i,j)~= 0
            if position_IadoptJ_in_rows(i,j)~=0
                DV_IJ_sparse = sparse(position_IadoptJ_in_rows(i,j),1,1,number_rows_for_DV(i,j),1);
            else
                DV_IJ_sparse = sparse(number_rows_for_DV(i,j),1);
            end
            DV_IJ_sparse_combine = [DV_IJ_sparse_combine; DV_IJ_sparse];
        else
        end
    end
end