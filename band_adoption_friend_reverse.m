clc
clear all

% adoptions = csvread('bandadoptions3.csv');               % note that adoptions has not subtracted 104
adoptions = csvread('bandadoptions_lenient_adopt.csv',1,0);    
adoptions = adoptions(:,1:3);
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
% friendlist = csvread('friends3.csv');
friendlist = csvread('new_friendlist_8088.csv',1,0);               % changes here !
bandtime = csvread('tbands3.csv');

bandtime(:,2) = bandtime(:,2)-old_bandt;
bandtime(:,3) = bandtime(:,3)-old_bandt;

timesplit(:,2) = timesplit(:,2)-old_bandt;
timesplit(:,3) = timesplit(:,3)-old_bandt;
timesplit(:,4) = timesplit(:,4)-old_bandt;

% largeNum_p = 68000000
largeNum_p = 120000000
% largeNum_p = 600000

member_p = zeros(largeNum_p,1);
friend_p = zeros(largeNum_p,1);
band_p = zeros(largeNum_p,1);
timeobs = zeros(largeNum_p,1);
DV = zeros(largeNum_p,1);
prob_adopt_week = zeros(largeNum_p,1);
% abs_week_diff = zeros(largeNum_p,1);
new_week_diff = zeros(largeNum_p,1);
A_week_ijt = zeros(largeNum_p,1);

ind = 1;
for i = 1:size(friendlist,1)                                         % i here  is the row num of friendlist
    member_id = friendlist(i,1);
    friend_id = friendlist(i,2);
    for j = 1:J
        adoptfij_time = band_adopt_mat(friend_id,j);
        adoptmij_time = band_adopt_mat(member_id,j);
%         if adoptfij_time~=0 && adoptmij_time~=0
        if adoptfij_time~=0 && adoptfij_time > timesplit(friend_id,3)
%             if adoptfij_time > timesplit(friend_id,3) && adoptmij_time > timesplit(member_id,3)
            if adoptmij_time > timesplit(member_id,3) || adoptmij_time==0
                pre_start = max([timesplit(member_id,3) bandtime(j,2)]);
                pre_end = min([timesplit(member_id,4) bandtime(j,3)]);        % still overlap period between friend obs and band obs ??
                interval = pre_end - pre_start +1;
                member_p(ind:ind+interval-1) = member_id;
                friend_p(ind:ind+interval-1) = friend_id;
                band_p(ind:ind+interval-1) = j;
                timeobs(ind:ind+interval-1) = pre_start:pre_end;
                if adoptmij_time~=0
                    DV(ind+(adoptmij_time-pre_start)) = 1;
                else
                end
%                 prob_adopt_week(ind:ind+interval-1) = prob_adoption_yearj(pre_start:pre_end,j);
                diffusion_remove_store = diffusion_jt(pre_start:pre_end,j)-DV(ind:ind+interval-1);      % remove adoption of the specific individual: cleaned bp
                prob_adopt_week(ind:ind+interval-1) = baseline_prob_smooth_func(diffusion_remove_store,10);
                week_diff = (pre_start:pre_end) - adoptfij_time;
                week_diff_rep = week_diff;
                %             abs_week_diff(ind:ind+interval-1) = abs(week_diff);
                week_diff(week_diff<0) = 0;
                new_week_diff(ind:ind+interval-1) = week_diff;
                week_diff_rep(week_diff_rep>=0) = 1;
                week_diff_rep(week_diff_rep<0) = 0;
                A_week_ijt(ind:ind+interval-1) = week_diff_rep;
                ind = ind+interval;
            else
            end
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
% abs_week_diff = abs_week_diff(1:(ind-1));
new_week_diff = new_week_diff(1:(ind-1));
A_week_ijt = A_week_ijt(1:(ind-1));

display('save as mat')

% matp = [member_p friend_p band_p timeobs DV prob_adopt_week abs_week_diff];
matp = [member_p friend_p band_p timeobs DV prob_adopt_week new_week_diff A_week_ijt];
% clearvars -EXCEPT matp
% save('matp2.mat','matp', '-v7.3') ;
save('matp_friend_reverse.mat','matp', '-v7.3') ;
% csvwrite('matp.csv',matp);

%% check number of band adoptions by members
load('matp_friend_reverse.mat')
mtind = find(matp(:,5));

TMAT = matp(1:63000,:);
mt_verif = zeros(length(mtind),1);
mtind = find(matp(:,5));
for i = 1:length(mtind)
    mt_temp = band_adopt_mat(matp(mtind(i),1),matp(mtind(i),3));
    mt_verif(i) = mt_temp==matp(mtind(i),4);
%     if i>4000
%         warning(msg)
%     end
end

