%% save for now
clc
clear

% load('matp_trend_rm.mat');
% load('innov_contin_trend_rm.mat');
% load('explor_contin_trend_rm.mat');
% load('matp_trend_lenient_rm.mat');
% load('innov_contin_trend_lenient_rm.mat');
% load('explor_contin_trend_lenient_rm.mat');
load('matp_trend_strict_rm.mat');
load('innov_contin_trend_strict_rm.mat');
load('explor_contin_trend_strict_rm.mat');

X = matp(:,1:6);                                    % trend on col 6
member_adopt_covariate = matp(:,7)>0;               
X = [X member_adopt_covariate];                     % member adopt on col 7
clearvars matp 

% load('matp_newbp_clean_former_rm.mat');
% load('matp_newbp_clean_lenient_rm.mat');
load('matp_newbp_clean_strict_rm.mat');
X = [X matp(:,6)];                                  % new bp on col 8
clearvars matp

matp = [X innov_contin explor_contin];
clearvars X innov_contin explor_contin

% ONLY for strict data
matp = matp(1:5091954,:);                   % quick & dirty temp bug fix for strict data
% clearvars cox_mat

%%
% load('matp_trend_rm.mat');

% last_obs = matp(1,3);
% band_change_ind = 1;
% for i= 2:size(matp,1)
%     current_obs = matp(i,3);
%     if current_obs~=last_obs
%         band_change_ind = [band_change_ind i];
%     else
%     end
%     last_obs = current_obs;
% end
% adopt_band_ids = matp(band_change_ind,3);             % with repeats   

last_obs = matp(1,2);
f_change_ind = 1;
for i= 2:size(matp,1)
    current_obs = matp(i,2);
    if current_obs~=last_obs
        f_change_ind = [f_change_ind i];
    else
    end
    last_obs = current_obs;
end
friend_ids = matp(f_change_ind,2);        

last_obs = matp(1,3);
start_ind = 1;
f_adopt = [];
adopt_band_ids = last_obs;
band_change_ind = 1;
for i= 2:size(matp,1)
    current_obs = matp(i,3);
    if current_obs~=last_obs || any(f_change_ind==i)
        f_adopt = [f_adopt sum(matp(start_ind:(i-1),5))];
        start_ind = i;
        adopt_band_ids = [adopt_band_ids current_obs];   % with repeats in band_ids
        band_change_ind = [band_change_ind i];
    else
    end
    last_obs = current_obs;
end
f_adopt = [f_adopt sum(matp(start_ind:i,5))];
f_adopt_rows = find(matp(:,5)==1);

cut_rows_after_adopt = [];
for j = 1:length(f_adopt_rows)
    cut_rows_after_adopt = [cut_rows_after_adopt band_change_ind(j):f_adopt_rows(j)];
end
save('cut_rows_after_adopt_strict.mat','cut_rows_after_adopt', '-v7.3') ;

length(f_adopt)                   % number of adoption observations in DV of the data matrix used for modeling
sum(f_adopt==1)
sum(matp(:,5))

sum(matp(band_change_ind,3) == adopt_band_ids')

adopt_ind = find(matp(:,5)==1);
sum(matp(adopt_ind,5))

TimeUntilEvent = adopt_ind - band_change_ind'+1;

%
cox_mat = zeros(size(matp,1),size(matp,2)+2);

rowsStore = zeros(length(band_change_ind)+1,1);
for i = 1:length(band_change_ind)-1
    observePeriod = matp(band_change_ind(i):band_change_ind(i+1),:);
    tStart = observePeriod(:,4)-observePeriod(1,4);
    tStop = tStart+1;
    temp = [observePeriod tStart tStop];
    untilEvent = temp(1:TimeUntilEvent(i),:);
    nrows = size(untilEvent,1);
    rowsStore(i+1) = rowsStore(i)+nrows;
    cox_mat(rowsStore(i)+1:rowsStore(i+1),:) = untilEvent;
end
observePeriod = matp(band_change_ind(i):end,:);
tStart = observePeriod(:,4)-observePeriod(1,4);
tStop = tStart+1;
temp = [observePeriod tStart tStop];
untilEvent = temp(1:TimeUntilEvent(i+1),:);
nrows = size(untilEvent,1);
rowsStore(i+2) = rowsStore(i+1)+nrows;
cox_mat(rowsStore(i+1)+1:rowsStore(i+2),:) = untilEvent;

cox_mat = cox_mat(1:rowsStore(end),:);

% save('cox_mat.mat','cox_mat', '-v7.3') ;
% csvwrite('cox_mat.csv',cox_mat);
% save('cox_mat_lenient.mat','cox_mat', '-v7.3') ;
% csvwrite('cox_mat_lenient.csv',cox_mat);
save('cox_mat_strict.mat','cox_mat', '-v7.3') ;
csvwrite('cox_mat_strict.csv',cox_mat);

% no week diff yet

%% test original dataset of old model used in thesis
load('matpstd2.mat');

temp = matp(matp(:,5)==1,:);
sum(temp(:,7)>0)                  % w: the original matp includes obs of friends adopt before members

last_obs = matp(1,3);
start_ind = 1;
f_adopt = [];
adopt_band_ids = last_obs;
band_change_ind = 1;
for i= 2:size(matp,1)
    current_obs = matp(i,3);
    if current_obs~=last_obs || any(f_change_ind==i)
        f_adopt = [f_adopt sum(matp(start_ind:(i-1),5))];
        start_ind = i;
        adopt_band_ids = [adopt_band_ids current_obs];   % with repeats in band_ids
        band_change_ind = [band_change_ind i];
    else
    end
    last_obs = current_obs;
end
f_adopt = [f_adopt sum(matp(start_ind:i,5))];

temp2 = matp(band_change_ind(2:end)-1,:);
sum(temp2(:,7)<=0)                 % w: the original matp also includes some bands adopted by friend but not by member, why?
