clc
clear

nm_6046_adopt = csvread('test_nm_full_artistID_adopt6046.csv',1,0);
fm_6046_adopt = csvread('test_fm_full_artistID_adopt6046.csv',1,0);

MAX_T = 527
J = max(fm_6046_adopt(:,2))                            % note this is >> 6046 !!
nm_I = max(nm_6046_adopt(:,1))
fm_I = max(fm_6046_adopt(:,1))
nm_adopt_mat = sparse(nm_6046_adopt(:,1),nm_6046_adopt(:,2),nm_6046_adopt(:,3),nm_I,J);
fm_adopt_mat = sparse(fm_6046_adopt(:,1),fm_6046_adopt(:,2),fm_6046_adopt(:,3),fm_I,J);

member_start_end = csvread('MEMBER_nm_fm_start_end_week.csv',1,0);
artist_start = csvread('ARTIST6046_start_week.csv',1,0);
NEW_AID_correspond_6046 = artist_start(:,1);

% get # of rows for DV
matrix_IJ_week_start = zeros(size(member_start_end,1),size(artist_start,1));
matrix_IJ_week_end = zeros(size(member_start_end,1),size(artist_start,1));
matrix_IadoptJ_week = zeros(size(member_start_end,1),size(artist_start,1));
for i = 1:size(member_start_end,1)
    for j = 1:size(artist_start,1)                % sequence(rows) of NEW_AID in artist_start correspond to 1-6046 band ids
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
sum(sum(matrix_IJ_week_end(find(matrix_IadoptJ_week>0))>=matrix_IadoptJ_week(find(matrix_IadoptJ_week>0))))
matrix_IJ_week_end(matrix_IadoptJ_week>0)=matrix_IadoptJ_week(matrix_IadoptJ_week>0);

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

% save('DV_IJ_sparse_combine.mat','DV_IJ_sparse_combine','-v7.3');        % save DV !!

load('DV_IJ_sparse_combine.mat')                     

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

% save('number_rows_for_DV.mat','number_rows_for_DV','-v7.3');
% save('start_index_rows_for_DV.mat','start_index_rows_for_DV','-v7.3');
% save('index_rows_for_DV.mat','index_rows_for_DV','-v7.3');
% save('matrix_IJ_week_start.mat','matrix_IJ_week_start','-v7.3');
% save('matrix_IJ_week_end.mat','matrix_IJ_week_end','-v7.3');
% save('matrix_IadoptJ_week.mat','matrix_IadoptJ_week','-v7.3');
% save('NEW_AID_correspond_6046.mat','NEW_AID_correspond_6046','-v7.3');

load('matrix_IadoptJ_week.mat');
load('number_rows_for_DV.mat');
load('start_index_rows_for_DV.mat');
load('index_rows_for_DV.mat');
load('matrix_IJ_week_start.mat');
load('matrix_IJ_week_end.mat');
load('NEW_AID_correspond_6046.mat');

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

% save('nm_adopt_cells.mat','nm_adopt_cells','-v7.3');
    
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

% save('fm_adopt_cells.mat','fm_adopt_cells','-v7.3');

load('nm_adopt_cells.mat')
load('fm_adopt_cells.mat');

% verify
very = find(~cellfun(@isempty,nm_adopt_cells));
sum(very == artist_start(:,1))
very = find(~cellfun(@isempty,fm_adopt_cells));
sum(very == artist_start(:,1))                             % verify that all non-empty cells (6046 in total) in nm_adopt_cells/fm_adopt_cells follows the exact order as in artist_start

nm_adopt_6046cells = nm_adopt_cells(~cellfun(@isempty, nm_adopt_cells));
fm_adopt_6046cells = fm_adopt_cells(~cellfun(@isempty, fm_adopt_cells));

% save('nm_adopt_6046cells.mat','nm_adopt_6046cells','-v7.3');
% save('fm_adopt_6046cells.mat','fm_adopt_6046cells','-v7.3');

load('nm_adopt_6046cells.mat')
load('fm_adopt_6046cells.mat');

% add similarity scores
load('Cosine_similarity_scores_neighbour_new_listen1_TFIDF.mat')
load('Cosine_user_uv_ind_neighbours_new_listens_TFIDF.mat');            % same as Cosine_user_uv_ind_neighbours_new_listen1_TFIDF
Cosine_user_uv_ind_nm = Cosine_user_uv_ind;
Cosine_similarity_score_nm = Cosine_similarity_scores;

load('Cosine_similarity_scores_friend_new_listen1_TFIDF.mat')
load('Cosine_user_uv_ind_friend_new_listen1_TFIDF.mat');
Cosine_user_uv_ind_fm = Cosine_user_uv_ind;
Cosine_similarity_score_fm = Cosine_similarity_scores;

% get the adoption times similarity with neighbours for each member M (M cells {J cells{all M's neighbours I*all time period Tij}})
% m_neighbour_cells = cell(size(member_start_end,1),1);
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

% % save('m_neighbour_cells.mat','m_neighbour_cells','-v7.3');
% save('similarity_m_neighbour_cells.mat','similarity_m_neighbour_cells','-v7.3');
% save('m_neighbour_adoption_multiply_similarity_cells.mat','m_neighbour_adoption_multiply_similarity_cells','-v7.3');
% save('m_neighbour_adoption_multiply_similarity_cells_Tij.mat','m_neighbour_adoption_multiply_similarity_cells_Tij','-v7.3');
% % clear m_neighbour_cells
% clear similarity_m_neighbour_cells
% clear m_neighbour_adoption_multiply_similarity_cells
% clear m_neighbour_adoption_multiply_similarity_cells_Tij


% load('m_neighbour_cells.mat')
load('similarity_m_neighbour_cells.mat');
% load('m_neighbour_adoption_multiply_similarity_cells.mat');
load('m_neighbour_adoption_multiply_similarity_cells_Tij.mat');

X_m_neighbour_adoption_times_similarity_combine = [];
for m = 1:size(matrix_IadoptJ_week,1)
    for j = 1:size(matrix_IadoptJ_week,2)
        if number_rows_for_DV(m,j)~= 0
            X_m_neighbour_adoption_times_similarity_sparse = m_neighbour_adoption_multiply_similarity_cells_Tij{m}{j}';
            X_m_neighbour_adoption_times_similarity_combine = ...
                [X_m_neighbour_adoption_times_similarity_combine; X_m_neighbour_adoption_times_similarity_sparse];
        else
        end
    end
end

% save('X_m_neighbour_adoption_times_similarity_combine.mat','X_m_neighbour_adoption_times_similarity_combine','-v7.3');
load('X_m_neighbour_adoption_times_similarity_combine.mat');

% get the adoption times similarity with friends for each member M (M cells {J cells{all M's friends I*all time period Tij}})
% m_neighbour_cells = cell(size(member_start_end,1),1);
similarity_m_friend_cells = cell(size(member_start_end,1),1);
m_friend_adoption_multiply_similarity_cells = cell(size(member_start_end,1),1);
m_friend_adoption_multiply_similarity_cellsD1 = cell(size(member_start_end,1),1);
m_friend_adoption_multiply_similarity_cellsD4 = cell(size(member_start_end,1),1);
m_friend_adoption_multiply_similarity_cells_Tij = cell(size(member_start_end,1),1);
m_friend_adoption_multiply_similarity_cells_TijD1 = cell(size(member_start_end,1),1);
m_friend_adoption_multiply_similarity_cells_TijD4 = cell(size(member_start_end,1),1);
for m = 1:size(member_start_end,1)
    m_friend_ID = Cosine_user_uv_ind_fm(find(Cosine_user_uv_ind_fm(:,1)==member_start_end(m,2)),2);
    similarity_m_friend_cells{m} = Cosine_similarity_score_fm(find(Cosine_user_uv_ind_fm(:,1)==member_start_end(m,2)));
%     m_neighbour_cells{m} = cell(size(artist_start,1),1);
    m_friend_adoption_multiply_similarity_cells{m} = cell(size(artist_start,1),1);
    m_friend_adoption_multiply_similarity_cellsD1{m} = cell(size(artist_start,1),1);
    m_friend_adoption_multiply_similarity_cellsD4{m} = cell(size(artist_start,1),1);
    m_friend_adoption_multiply_similarity_cells_Tij{m} = cell(size(artist_start,1),1);
    m_friend_adoption_multiply_similarity_cells_TijD1{m} = cell(size(artist_start,1),1);
    m_friend_adoption_multiply_similarity_cells_TijD4{m} = cell(size(artist_start,1),1);
    for j = 1:size(artist_start,1)
%         m_neighbour_cells{m}{j} = nm_adopt_6046cells{j}(m_neighbours_ID,:);
        m_friend_cells_mj = fm_adopt_6046cells{j}(m_friend_ID,:);
        temp = similarity_m_friend_cells{m}'*m_friend_cells_mj;
        m_friend_adoption_multiply_similarity_cells{m}{j} = sparse(1,find(temp),temp(find(temp)),1,MAX_T);
        m_friend_adoption_multiply_similarity_cellsD1{m}{j} = [sparse(0,[],[],1,1) m_friend_adoption_multiply_similarity_cells{m}{j}(1:end-1)];         % adopted 1 week ago
        m_friend_adoption_multiply_similarity_cellsD4{m}{j} = [sparse(0,[],[],1,4) m_friend_adoption_multiply_similarity_cells{m}{j}(1:end-4)];         % adopted 4 week ago
        if matrix_IJ_week_start(m,j)>matrix_IJ_week_end(m,j)              % for the bands appear after member becomes inactive, do not include
            m_friend_adoption_multiply_similarity_cells_Tij{m}{j} = [];
            m_friend_adoption_multiply_similarity_cells_TijD1{m}{j} = [];
            m_friend_adoption_multiply_similarity_cells_TijD4{m}{j} = [];
        elseif matrix_IadoptJ_week(m,j)>0                                 % if not truncated
            m_friend_adoption_multiply_similarity_cells_Tij{m}{j} = ...
                m_friend_adoption_multiply_similarity_cells{m}{j}(matrix_IJ_week_start(m,j):min(matrix_IJ_week_end(m,j),matrix_IadoptJ_week(m,j)));
            m_friend_adoption_multiply_similarity_cells_TijD1{m}{j} = ...
                m_friend_adoption_multiply_similarity_cellsD1{m}{j}(matrix_IJ_week_start(m,j):min(matrix_IJ_week_end(m,j),matrix_IadoptJ_week(m,j)));
            m_friend_adoption_multiply_similarity_cells_TijD4{m}{j} = ...
                m_friend_adoption_multiply_similarity_cellsD4{m}{j}(matrix_IJ_week_start(m,j):min(matrix_IJ_week_end(m,j),matrix_IadoptJ_week(m,j)));
        else
            m_friend_adoption_multiply_similarity_cells_Tij{m}{j} = ...
                m_friend_adoption_multiply_similarity_cells{m}{j}(matrix_IJ_week_start(m,j):matrix_IJ_week_end(m,j));
            m_friend_adoption_multiply_similarity_cells_TijD1{m}{j} = ...
                m_friend_adoption_multiply_similarity_cellsD1{m}{j}(matrix_IJ_week_start(m,j):matrix_IJ_week_end(m,j));
            m_friend_adoption_multiply_similarity_cells_TijD4{m}{j} = ...
                m_friend_adoption_multiply_similarity_cellsD4{m}{j}(matrix_IJ_week_start(m,j):matrix_IJ_week_end(m,j));
        end
    end
end

% save('similarity_m_friend_cells.mat','similarity_m_friend_cells','-v7.3');
% save('m_friend_adoption_multiply_similarity_cells.mat','m_friend_adoption_multiply_similarity_cells','-v7.3');
% save('m_friend_adoption_multiply_similarity_cells_Tij.mat','m_friend_adoption_multiply_similarity_cells_Tij','-v7.3');
% save('m_friend_adoption_multiply_similarity_cellsD1.mat','m_friend_adoption_multiply_similarity_cellsD1','-v7.3');
% save('m_friend_adoption_multiply_similarity_cells_TijD1.mat','m_friend_adoption_multiply_similarity_cells_TijD1','-v7.3');
% save('m_friend_adoption_multiply_similarity_cellsD4.mat','m_friend_adoption_multiply_similarity_cellsD4','-v7.3');
% save('m_friend_adoption_multiply_similarity_cells_TijD4.mat','m_friend_adoption_multiply_similarity_cells_TijD4','-v7.3');

clear similarity_m_friend_cells
clear m_friend_adoption_multiply_similarity_cells
clear m_friend_adoption_multiply_similarity_cells_Tij
clear m_friend_adoption_multiply_similarity_cellsD1
clear m_friend_adoption_multiply_similarity_cells_TijD1
clear m_friend_adoption_multiply_similarity_cellsD4
clear m_friend_adoption_multiply_similarity_cells_TijD4

load('similarity_m_friend_cells.mat');
load('m_friend_adoption_multiply_similarity_cells_Tij.mat');
load('m_friend_adoption_multiply_similarity_cells_TijD1.mat');
load('m_friend_adoption_multiply_similarity_cells_TijD4.mat');

X_m_friend_adoption_times_similarity_combine = [];
X_m_friend_adoption_times_similarity_combineD1 = [];
X_m_friend_adoption_times_similarity_combineD4 = [];
for m = 1:size(matrix_IadoptJ_week,1)
    for j = 1:size(matrix_IadoptJ_week,2)
        if number_rows_for_DV(m,j)~= 0
            X_m_friend_adoption_times_similarity_sparse = m_friend_adoption_multiply_similarity_cells_Tij{m}{j}';
            X_m_friend_adoption_times_similarity_sparseD1 = m_friend_adoption_multiply_similarity_cells_TijD1{m}{j}';
            X_m_friend_adoption_times_similarity_sparseD4 = m_friend_adoption_multiply_similarity_cells_TijD4{m}{j}';
            X_m_friend_adoption_times_similarity_combine = ...
                [X_m_friend_adoption_times_similarity_combine; X_m_friend_adoption_times_similarity_sparse];
            X_m_friend_adoption_times_similarity_combineD1 = ...
                [X_m_friend_adoption_times_similarity_combineD1; X_m_friend_adoption_times_similarity_sparseD1];
            X_m_friend_adoption_times_similarity_combineD4 = ...
                [X_m_friend_adoption_times_similarity_combineD4; X_m_friend_adoption_times_similarity_sparseD4];
        else
        end
    end
end

% save('X_m_friend_adoption_times_similarity_combine.mat','X_m_friend_adoption_times_similarity_combine','-v7.3');
% save('X_m_friend_adoption_times_similarity_combineD1.mat','X_m_friend_adoption_times_similarity_combineD1','-v7.3');
% save('X_m_friend_adoption_times_similarity_combineD4.mat','X_m_friend_adoption_times_similarity_combineD4','-v7.3');
load('X_m_friend_adoption_times_similarity_combine.mat');
load('X_m_friend_adoption_times_similarity_combineD1.mat');
load('X_m_friend_adoption_times_similarity_combineD4.mat');

%%
% include baseline probabilities and google trends variables
% baseline probabilities computed based on all listens and adoptions from lastfm dataset 
baselineProbListens = csvread('dt7_member_plus_friend_neighbour_names_fix_weeks_artist6046ID_union_week_listens.csv',1,0);
baselineProbAdopts = csvread('test_fm_nm_artistID6046_adopt_week_baseline.csv',1,0);

% mean(baselineProbAdopts(:,3))
% mean(baselineProbListens(:,3))
                   
baselineProbListens_mat = sparse(baselineProbListens(:,1),baselineProbListens(:,2),baselineProbListens(:,3),J,MAX_T);
baselineProbAdopts_mat = sparse(baselineProbAdopts(:,1),baselineProbAdopts(:,2),baselineProbAdopts(:,3),J,MAX_T);

X_baselineProbListens = zeros(size(DV_IJ_sparse_combine));
X_baselineProbAdopts = zeros(size(DV_IJ_sparse_combine));

for m = 1:size(matrix_IadoptJ_week,1)
    for j = 1:size(matrix_IadoptJ_week,2)
        if number_rows_for_DV(m,j)~= 0
            X_baselineProbListens(start_index_rows_for_DV(m,j):index_rows_for_DV(m,j)) = ...
                baselineProbListens_mat(NEW_AID_correspond_6046(j), matrix_IJ_week_start(m,j):matrix_IJ_week_end(m,j));
            X_baselineProbAdopts(start_index_rows_for_DV(m,j):index_rows_for_DV(m,j)) = ...
                baselineProbAdopts_mat(NEW_AID_correspond_6046(j), matrix_IJ_week_start(m,j):matrix_IJ_week_end(m,j));
        else
        end
    end
end

% save('X_baselineProbListens.mat','X_baselineProbListens','-v7.3');
% save('X_baselineProbAdopts.mat','X_baselineProbAdopts','-v7.3');

load('X_baselineProbListens.mat');
load('X_baselineProbAdopts.mat');

% view correlations between baseline probabilities computed from adoptions vs listens 
corr(X_baselineProbAdopts, X_baselineProbListens)

% verify that the BAND_IDs used to collect google trends for the 6046 bands ...
% has the same order as the sequence of NEW_AIDs used to construct the 165*6046 matrix 
old_new_band_id_trans = csvread('artist_new_old_id_6046.csv',1,0);
sum(old_new_band_id_trans(:,2)==NEW_AID_correspond_6046)

% include google trends collected before
predict_trend_log4061_original = csvread('predict_trend_log4061_original.csv',1,0);
predict_trend_log4061_lenient = csvread('predict_trend_log4061_lenient.csv',1,0);

% min(min(matrix_IJ_week_start))
% max(predict_trend_log4061_original(:,2))
% Our new bands data in this script are coded between 102 to 527 weeks, but the
% google trends data are coded between 1-423 weeks, thus need to address google 
% trends data by adding 104 weeks
predict_trend_log4061_original(:,2) = predict_trend_log4061_original(:,2)+104;
predict_trend_log4061_lenient(:,2) = predict_trend_log4061_lenient(:,2)+104;

bandlist4061 = unique(predict_trend_log4061_lenient(:,1));
% min(abs(predict_trend_log4061_original(find(predict_trend_log4061_original(:,3)~=0),3)))
predict_trend_log4061_original(find(predict_trend_log4061_original(:,3)==0),3)=0.00000001;

predict_trend_log4061_original_mat = sparse(predict_trend_log4061_original(:,1),predict_trend_log4061_original(:,2),predict_trend_log4061_original(:,3),6046,MAX_T);
predict_trend_log4061_lenient_mat = sparse(predict_trend_log4061_lenient(:,1),predict_trend_log4061_lenient(:,2),predict_trend_log4061_lenient(:,3),6046,MAX_T);

X_predict_trend_log4061_original = zeros(size(DV_IJ_sparse_combine));
X_predict_trend_log4061_lenient = zeros(size(DV_IJ_sparse_combine));

for m = 1:size(matrix_IadoptJ_week,1)
    for j = 1:size(matrix_IadoptJ_week,2)
        if number_rows_for_DV(m,j)~= 0
            X_predict_trend_log4061_original(start_index_rows_for_DV(m,j):index_rows_for_DV(m,j)) = ...
                predict_trend_log4061_original_mat(j, matrix_IJ_week_start(m,j):matrix_IJ_week_end(m,j));
            X_predict_trend_log4061_lenient(start_index_rows_for_DV(m,j):index_rows_for_DV(m,j)) = ...
                predict_trend_log4061_lenient_mat(j, matrix_IJ_week_start(m,j):matrix_IJ_week_end(m,j));
        else
        end
    end
end

corr(X_predict_trend_log4061_original, X_baselineProbListens)
corr(X_predict_trend_log4061_lenient, X_baselineProbListens)
corr(X_predict_trend_log4061_original, X_baselineProbAdopts)
corr(X_predict_trend_log4061_lenient, X_baselineProbAdopts)

corr(X_predict_trend_log4061_original(X_predict_trend_log4061_original>0), X_baselineProbListens(X_predict_trend_log4061_original>0))
corr(X_predict_trend_log4061_lenient(X_predict_trend_log4061_lenient>0), X_baselineProbListens(X_predict_trend_log4061_lenient>0))
corr(X_predict_trend_log4061_original(X_predict_trend_log4061_original>0), X_baselineProbAdopts(X_predict_trend_log4061_original>0))
corr(X_predict_trend_log4061_lenient(X_predict_trend_log4061_lenient>0), X_baselineProbAdopts(X_predict_trend_log4061_lenient>0))

% X_predict_trend_log4061_original has higher corr with both X_baselineProbListens and X_baselineProbListens

% save('X_predict_trend_log4061_original.mat','X_predict_trend_log4061_original','-v7.3');
% save('X_predict_trend_log4061_lenient.mat','X_predict_trend_log4061_lenient','-v7.3');

load('X_predict_trend_log4061_original.mat')
load('X_predict_trend_log4061_lenient.mat')
