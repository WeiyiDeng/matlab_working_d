clc
clear

nm_6046_adopt = csvread('test_nm_full_artistID_adopt6046.csv',1,0);
fm_6046_adopt = csvread('test_fm_full_artistID_adopt6046.csv',1,0);

J = max(fm_6046_adopt(:,2))                            % note this is >> 6046 !!
nm_I = max(nm_6046_adopt(:,1))
fm_I = max(fm_6046_adopt(:,1))
nm_adopt_mat = sparse(nm_6046_adopt(:,1),nm_6046_adopt(:,2),nm_6046_adopt(:,3),nm_I,J);
fm_adopt_mat = sparse(fm_6046_adopt(:,1),fm_6046_adopt(:,2),fm_6046_adopt(:,3),fm_I,J);

member_start_end = csvread('MEMBER_nm_fm_start_end_week.csv',1,0);
artist_start = csvread('ARTIST6046_start_week.csv',1,0);

% get # of rows for DV
matrix_IJ_week_start = zeros(size(member_start_end,1),size(artist_start,1));
matrix_IJ_week_end = zeros(size(member_start_end,1),size(artist_start,1));
matrix_IadoptJ_week = zeros(size(member_start_end,1),size(artist_start,1));
for i = 1:size(member_start_end,1)
    for j = 1:size(artist_start,1)
        matrix_IJ_week_start(i,j) = max(artist_start(j,2),member_start_end(i,3));
        matrix_IJ_week_end(i,j) = member_start_end(i,4);
        member_ID_nm = member_start_end(i,1);
        artist_ID = artist_start(j,1);
        matrix_IadoptJ_week(i,j) = nm_adopt_mat(member_ID_nm,artist_ID);
    end
end

% verify
sum(sum(matrix_IJ_week_end<matrix_IJ_week_start))      % some bands appear after member becomes inactive

% sum(sum(matrix_IJ_week_end<matrix_IadoptJ_week))

% verify
sum(sum(matrix_IJ_week_start<=matrix_IadoptJ_week)) == sum(sum(matrix_IadoptJ_week>0))     % verify # of adoptions in DV

number_rows_for_DV = matrix_IJ_week_end-matrix_IJ_week_start+1;
number_rows_for_DV(number_rows_for_DV<0) = 0;          % set bands appear after member becomes inactive to have 0 rows

position_IadoptJ_in_rows = matrix_IadoptJ_week-matrix_IJ_week_start+1;
sum(sum(position_IadoptJ_in_rows>0))                   % verify once more
position_IadoptJ_in_rows(position_IadoptJ_in_rows<0) = 0;

% verify
sum(sum(position_IadoptJ_in_rows+matrix_IJ_week_start-1>matrix_IJ_week_end))
sum(sum((matrix_IJ_week_end<matrix_IJ_week_start).*...
    (position_IadoptJ_in_rows+matrix_IJ_week_start-1>matrix_IJ_week_end)))       % verfiy all these bands do not add right belong to the bands that appear after member becomes inactive

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

save('DV_IJ_sparse_combine.mat','DV_IJ_sparse_combine','-v7.3');

% verify # of adoptions in DV
sum(DV_IJ_sparse_combine)

trans_number_rows_for_DV = number_rows_for_DV';
index_rows_for_DV = reshape(cumsum(trans_number_rows_for_DV(:)),size(number_rows_for_DV,2),size(number_rows_for_DV,1));
index_rows_for_DV = index_rows_for_DV';

temp_trans = index_rows_for_DV';
temp = temp_trans(:);
temp = [0; temp];
temp(end) = [];
temp = temp+1;
start_index_rows_for_DV = reshape(temp,size(number_rows_for_DV,2),size(number_rows_for_DV,1));
start_index_rows_for_DV = start_index_rows_for_DV';

save('number_rows_for_DV.mat','number_rows_for_DV','-v7.3');
save('start_index_rows_for_DV.mat','start_index_rows_for_DV','-v7.3');
save('index_rows_for_DV.mat','index_rows_for_DV','-v7.3');

% DV construction finished !

% construct adoption cells
nm_adopt_cells = cell(max(nm_6046_adopt(:,2)),1);
for j = 1:max(nm_6046_adopt(:,2))
% for j = 1:109
    temp = find(nm_6046_adopt(:,2)==j);
    if isempty(temp)
        nm_adopt_cells{j} = [];
    else
        j_matrix = nm_6046_adopt(temp,:);
        nm_adopt_cells{j} = sparse(j_matrix(:,1),j_matrix(:,3),1,max(nm_6046_adopt(:,1)),max(nm_6046_adopt(:,3)));
    end
end

save('nm_adopt_cells.mat','nm_adopt_cells','-v7.3');
    
fm_adopt_cells = cell(max(fm_6046_adopt(:,2)),1);
for j = 1:max(fm_6046_adopt(:,2))
% for j = 1:109
    temp = find(fm_6046_adopt(:,2)==j);
    if isempty(temp)
        fm_adopt_cells{j} = [];
    else
        j_matrix = nm_6046_adopt(temp,:);
        fm_adopt_cells{j} = sparse(j_matrix(:,1),j_matrix(:,3),1,max(fm_6046_adopt(:,1)),max(fm_6046_adopt(:,3)));
    end
end

save('fm_adopt_cells.mat','fm_adopt_cells','-v7.3');

% verify
very = find(~cellfun(@isempty,nm_adopt_cells));
sum(very == artist_start(:,1))
very = find(~cellfun(@isempty,fm_adopt_cells));
sum(very == artist_start(:,1))                             % verify that all non-empty cells (6046 in total) in nm_adopt_cells/fm_adopt_cells follows the exact order as in artist_start

nm_adopt_6046cells = nm_adopt_cells(~cellfun(@isempty, nm_adopt_cells));
fm_adopt_6046cells = fm_adopt_cells(~cellfun(@isempty, fm_adopt_cells));

save('nm_adopt_6046cells.mat','nm_adopt_6046cells','-v7.3');
save('fm_adopt_6046cells.mat','fm_adopt_6046cells','-v7.3');
