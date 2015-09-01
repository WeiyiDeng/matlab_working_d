clc
clear all

adoptions = csvread('bandadoptions3.csv');               % note that adoptions has not subtracted 104
adopt_ones = ones(size(adoptions, 1),1);

% total number of weeks: length(105:527) = 423
T = 423
I = 8320                                                 % 165 members
J = 6046
old_bandt = 104
band_adoption = cell(T,1);

for t = 1:T
    [r c v] = find(adoptions(:,3)==(t+old_bandt));    
    adopt_ind = adoptions(r,:);
    band_adoption{t} = sparse(adopt_ind(:,1),adopt_ind(:,2),adopt_ones(r),I,J);
end

band_adopt_mat = sparse(adoptions(:,1),adoptions(:,2),adoptions(:,3)-old_bandt,I,J);                % important fix !!

% wwwwwwwwwwwwwwwwwwwwww!!!!!!!!!!!!!!!!!!!!!!!!!!!!


% % test
% sum_a = 0;
% for t = 1:T
%     sum_a = sum_a + sum(sum([band_adoption{t}]));
% end
% sum_a

diffusion_jt = zeros(T,J);
for t = 1:T
    diffusion_jt(t,:) = sum(band_adoption{t},1);           % T*J
end

% overall probability of adopting band j at time t (pdf)
sum_adoptions_j = sum(diffusion_jt,1);
% prob_adoption_tj = diffusion_jt./repmat(sum_adoptions_j,T,1);    % PAjt 

diffusion_year = zeros(T,J);
for i = 1:7
    ind = (i-1)*52+1;
    year_adopt = sum(diffusion_jt(ind:ind+51,:),1);
    diffusion_year(ind:ind+51,:) = repmat(year_adopt,52,1);
end
year_adopt_8th = sum(diffusion_jt(365:423,:),1);
diffusion_year(365:423,:) = repmat(year_adopt_8th,59,1);

% sum_adoptions_yearj = sum(diffusion_year,1);
prob_adoption_yearj = diffusion_year./repmat(sum_adoptions_j,T,1);    % PAjt 

cumu_diffusion = cumsum(diffusion_jt);                    % T*J
 
timesplit = csvread('tsplit3.csv');
friendlist = csvread('friends3.csv');
bandtime = csvread('tbands3.csv');

bandtime(:,2) = bandtime(:,2)-old_bandt;
bandtime(:,3) = bandtime(:,3)-old_bandt;

timesplit(:,2) = timesplit(:,2)-old_bandt;
timesplit(:,3) = timesplit(:,3)-old_bandt;
timesplit(:,4) = timesplit(:,4)-old_bandt;

largeNum_p = 68000000
% largeNum_p = 600000

member_p = zeros(largeNum_p,1);
friend_p = zeros(largeNum_p,1);
band_p = zeros(largeNum_p,1);
timeobs = zeros(largeNum_p,1);
DV = zeros(largeNum_p,1);
prob_adopt_week = zeros(largeNum_p,1);
abs_week_diff = zeros(largeNum_p,1);


ind = 1;
for i = 1:size(friendlist,1)                                         % i here  is the row num of friendlist
    member_id = friendlist(i,1);
    friend_id = friendlist(i,2);
    for j = 1:J
        adoptfij_time = band_adopt_mat(friend_id,j);
        adoptmij_time = band_adopt_mat(member_id, j);
        if adoptfij_time~=0 && adoptmij_time~=0
            pre_start = max([timesplit(friend_id,2) bandtime(j,2)]);
            pre_end = min([timesplit(friend_id,4) bandtime(j,3)]);
            interval = pre_end - pre_start +1;
            member_p(ind:ind+interval-1) = member_id;
            friend_p(ind:ind+interval-1) = friend_id;
            band_p(ind:ind+interval-1) = j;
            timeobs(ind:ind+interval-1) = pre_start:pre_end;
            DV(ind+(adoptfij_time-pre_start)) = 1;
            prob_adopt_week(ind:ind+interval-1) = prob_adoption_yearj(pre_start:pre_end,j);
            abs_week_diff(ind:ind+interval-1) = abs((pre_start:pre_end) - adoptmij_time);
            ind = ind+interval;
        else
        end
    end
%     if i >=55
%         break
%     end
end

member_p = member_p(1:(ind-1));
friend_p = friend_p(1:(ind-1));
band_p = band_p(1:(ind-1));
timeobs = timeobs(1:(ind-1));
DV = DV(1:(ind-1));
prob_adopt_week = prob_adopt_week(1:(ind-1));
abs_week_diff = abs_week_diff(1:(ind-1));

display('save as mat')

matp = [member_p friend_p band_p timeobs DV prob_adopt_week abs_week_diff];
% clearvars -EXCEPT matp
save('matp.mat','matp', '-v7.3') ;
% csvwrite('matp.csv',matp);

%% member rows
row_ind = 1
row_mid = [];
row_num = [];
new_row = matp(1,1);
for i = 2:size(matp,1)
    old_row = new_row;
    new_row = matp(i,1);
    if new_row == old_row
        row_ind = row_ind + 1;
    else
        row_mid = [row_mid old_row];
        row_num = [row_num row_ind];
        row_ind = 1;
    end
end
row_mid = [row_mid new_row];
row_num = [row_num row_ind];

sum(row_num)

save('row_num.mat','row_num');
save('row_mid.mat','row_mid');
% csvwrite('row_num.csv',row_num);
% csvwrite('row_mid.csv',row_mid);

%% create member dummies
load('matp.mat');
matp_len = size(matp,1);
week_diff_for_d = matp(:,7);
clearvars matp

load('row_num.mat')
load('row_mid.mat')
members_for_d = csvread('members_for_dummies.csv');

row_cumsum = cumsum(row_num,2);

% w = sparse([3:5 7],ones(1,4).*2,1,10,5);
% ww = repmat(w,2,1);                % after repmat still a sparse matrix
% www = w*rand(5);                   % after matrix mutiplication becomes a dense matrix 

members_rind = [];
dummies_cind = [];
num_ind = 1
for i = 1:size(members_for_d,1)
    [r c v] = find(row_mid == members_for_d(i,1));
    start_ind = row_cumsum(c-1)+1;
    end_ind = row_cumsum(c);
    members_rind = [members_rind start_ind:end_ind];
    dummies_cind = [dummies_cind ones(1,length(start_ind:end_ind)).*num_ind];
    num_ind = num_ind+1;
end

member_dummies = sparse(members_rind,dummies_cind,1,matp_len,size(members_for_d,1));

save('member_dummies.mat','member_dummies');

member_dummies_week_d = member_dummies;
for i = 1:size(member_dummies_week_d,2)
    member_dummies_week_d(:,i) = member_dummies_week_d(:,i).*week_diff_for_d;
end

save('member_dummies_week_d.mat','member_dummies_week_d');

%% friend and band obs overlaps
row_ind = 1
row_fid = [];
row_interval = [];
new_row = matp(1,1);
for i = 2:size(matp,1)
    old_row = new_row;
    new_row = matp(i,2);
    if new_row == old_row
        row_ind = row_ind + 1;
    else
        row_fid = [row_fid old_row];
        row_interval = [row_interval row_ind];
        row_ind = 1;
    end
end
row_fid = [row_fid new_row];
row_interval = [row_interval row_ind];

save('row_fid.mat','row_fid');
save('row_interval.mat','row_interval');