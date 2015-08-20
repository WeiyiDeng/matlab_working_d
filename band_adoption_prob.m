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

cumu_diffusion = cumsum(diffusion_jt);                    % T*J

% overall probability of adopting band j at time t (pdf)
sum_adoptions_j = sum(diffusion_jt,1);
prob_adoption_tj = diffusion_jt./repmat(sum_adoptions_j,T,1);    % PAjt   

timesplit = csvread('tsplit3.csv');
friendlist = csvread('friends3.csv');
bandtime = csvread('tbands3.csv');

bandtime(:,2) = bandtime(:,2)-old_bandt;
bandtime(:,3) = bandtime(:,3)-old_bandt;

timesplit(:,2) = timesplit(:,2)-old_bandt;
timesplit(:,3) = timesplit(:,3)-old_bandt;
timesplit(:,4) = timesplit(:,4)-old_bandt;

% largeNum_p = 80000000
largeNum_p = 600000

member_p = zeros(largeNum_p,1);
friend_p = zeros(largeNum_p,1);
band_p = zeros(largeNum_p,1);
timeobs = zeros(largeNum_p,1);
DV = zeros(largeNum_p,1);
prob_adopt_week = zeros(largeNum_p,1);
abs_week_diff = zeros(largeNum_p,1);

% rownum = zeros(size(friendlist,1),2);
% row_ind = 1;

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
            prob_adopt_week(ind:ind+interval-1) = prob_adoption_tj(pre_start:pre_end,j);
            abs_week_diff(ind:ind+interval-1) = abs((pre_start:pre_end) - adoptmij_time);
            ind = ind+interval;
        else
        end
    end
%     rownum(row_ind,1) = member_id;
%     rownum(row_ind,2) = ind-1;
%     row_ind = row_ind +1;
    if i >=55
        break
    end
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








