clc
clear

% load('cox_mat_strict.mat');
% load('cox_mat_lenient.mat');
load('cox_mat.mat');

userID_1st_listen = csvread('userID_1st_listen.csv');       % total number of weeks 527
bandID_1st_listen = csvread('bandID_1st_listen.csv');       % total number of weeks 527

userID_1st_listen(:,2) = userID_1st_listen(:,2)-104;
bandID_1st_listen(:,2) = bandID_1st_listen(:,2)-104;

user_first_listen = zeros(size(cox_mat,1),1);
band_first_listen = zeros(size(cox_mat,1),1);

for i = 1:length(user_first_listen)
    user_first_listen(i) = userID_1st_listen(cox_mat(i,2),2);
end

for j = 1:length(band_first_listen)
    band_first_listen(j) = bandID_1st_listen(cox_mat(j,3),2);
end

user_band_first_listen_smallest = user_first_listen.*(user_first_listen > band_first_listen)...        % later one
    +band_first_listen.*(user_first_listen <= band_first_listen);

sum(cox_mat(:,4)<user_band_first_listen_smallest)              % checking

diff_weeks = cox_mat(:,4)-user_band_first_listen_smallest;

last_obs = cox_mat(1,3);
last_friend = cox_mat(1,2);
band_change_ind = 1;
for i= 2:size(cox_mat,1)
    current_obs = cox_mat(i,3);
    current_friend = cox_mat(i,2);
    if current_obs~=last_obs
        band_change_ind = [band_change_ind i];
    else
    end
    last_obs = current_obs;
end
band_ids = cox_mat(band_change_ind,3);   
obs_original_start_week = cox_mat(band_change_ind,4);

obs_new_start_week = user_band_first_listen_smallest(band_change_ind);

diff_old_new_start_week = obs_original_start_week-obs_new_start_week;

tStart_changed = zeros(size(cox_mat,1),1);
tEnd_changed = zeros(size(cox_mat,1),1);
for i = 1:(length(band_change_ind)-1)
    period_to_change = band_change_ind(i):(band_change_ind(i+1)-1);
    tStart_changed(period_to_change) = cox_mat(period_to_change,21)+diff_old_new_start_week(i);
    tEnd_changed(period_to_change) = cox_mat(period_to_change,22)+diff_old_new_start_week(i);
end
period_to_change = band_change_ind(i+1):size(cox_mat,1);
tStart_changed(period_to_change) = cox_mat(period_to_change,21)+diff_old_new_start_week(i+1);
tEnd_changed(period_to_change) = cox_mat(period_to_change,22)+diff_old_new_start_week(i+1);

max(tStart_changed)                               % checking
max(tEnd_changed)
min(tStart_changed)
min(tEnd_changed)

cox_mat_t = [cox_mat(:,1:20) tStart_changed tEnd_changed];

% save('cox_mat_strict_tStart.mat','cox_mat_t','-v7.3');
% csvwrite('cox_mat_strict_tStart.csv',cox_mat_t);

% save('cox_mat_lenient_tStart.mat','cox_mat_t','-v7.3');
% csvwrite('cox_mat_lenient_tStart.csv',cox_mat_t);

save('cox_mat_tStart.mat','cox_mat_t','-v7.3');
csvwrite('cox_mat_tStart.csv',cox_mat_t);

% review data to see whether data definition criteria results in systematic structure in the distribution of T 
load('cox_mat_tStart.mat');
start_T_ind = cox_mat_t(band_change_ind, 21);
hist(start_T_ind)
mean(start_T_ind)
min(start_T_ind)
max(start_T_ind)

