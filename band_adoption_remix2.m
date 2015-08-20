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

band_pop = zeros(T,J);
band_pop(1:8,:) = cumu_diffusion(1:8,:);
for t = 9:T                             % 8+1
    band_pop(t,:) = cumu_diffusion(t,:)- cumu_diffusion(t-8,:);
end


% test: 
% max(max(cumu_diffusion))
% [row,col]=find(cumu_diffusion==108)     % MGMT BAND_ID=4096
% cumu_diffusion(:,150)                   % Adele
Timew = 1:T;
scatter(Timew,cumu_diffusion(:,4096));
scatter(Timew,cumu_diffusion(:,150));
%%
timesplit = csvread('tsplit3.csv');

% % test: these people should already be those with new band listens
% find(timesplit(:,2:4)<=old_bandt);

% explorer = zeros(I,1);
% for i = 1:I
%     for t = old_bandt+1:timesplit(i,3)             
%         explorer(i) = explorer(i)+sum(band_adoption{t-old_bandt}(i,:));
%     end
% end
% 
% csvwrite('explorer3.csv',explorer);

explorer = csvread('explorer3.csv');

% test: moving average of 8 weeks
% Timeww = 1:T-8;
% for t = 1:T-8
%     wtry(t) = mean(diffusion_jt(t:t+7,150));
% end
% plot(wtry)
% scatter(Timeww,wtry)
% pk = find(wtry==max(wtry))
% pkw = pk(1)+3

% 
% moving_avg = 8;
% peakj = zeros(J,1);
% for j = 1:J
%     moving_w = zeros(1,T-moving_avg);
%     for t = 1:T-moving_avg
%         moving_w(t) = mean(diffusion_jt(t:t+moving_avg-1,j));
%     end
%     peakw = find(moving_w==max(moving_w));
%     if length(peakw) ==1
%     peakj(j) = peakw + moving_avg/2-1;     % there are still multiple max peak adoption weeks, how to fix?
%     else
%         ind_pkw = ceil(length(peakw)/2);
%         peakj(j) = peakw(ind_pkw) + moving_avg/2;         %-1;
%     end
% end

% introdate = csvread('introdatej3.csv');
% 
% introdate(:,2) = introdate(:,2)-old_bandt;    

% % test: some have same peak week and introweek
% % E: J = 139, use diffusion_jt(:,139) to see week adoptions
% % introdate = peakw = 241
% testA = peakj-introdate(:,2);
% find(testA<0);

% peakj = [introdate(:,1) peakj];
% adoptions(:,3) = adoptions(:,3)-old_bandt;
% csvwrite('peak3.csv',peakj);
% csvwrite('introdate3.csv',introdate);
% csvwrite('mod_adoptions3.csv',adoptions);
% give headers "BAND_ID", "USER_ID", "week_mod" etc.
% to the csv files in Emeditor before importing to SQL 

% adoptingweek = csvread('adoptingweek3.csv');

% EAij = 1-(adoptingweek(:,3)-adoptingweek(:,4))./(adoptingweek(:,5)-adoptingweek(:,4));
% % some have the same intro week and peak week, get NAN for 0 in denominator
% % test: 
% testC = adoptingweek(:,5)==adoptingweek(:,4);
% find(testC);

% ind_beforepk = EAij > 0;                  % also get rid of NANs here
% ind_beforepk = logical(1-ind_beforepk);
% EAij(ind_beforepk) = 0;
% 
% EAi = zeros(I,1);
% for i = 1:I
%     ind_i = find(adoptingweek(:,1)==i);
%     EAi(i) = mean(EAij(ind_i));           % changed from median to mean !!
% end
% 
% EAi = [(1:I)' EAi];
% csvwrite('EAi3.csv',EAi);

EAcat = csvread('EAdummy3.csv');          % dummy created from excel using EAcats.csv

%%
% friend trials
friendlist = csvread('friends3.csv');
bandtime = csvread('tbands3.csv');

bandtime(:,2) = bandtime(:,2)-old_bandt;
bandtime(:,3) = bandtime(:,3)-old_bandt;

timesplit(:,2) = timesplit(:,2)-old_bandt;
timesplit(:,3) = timesplit(:,3)-old_bandt;
timesplit(:,4) = timesplit(:,4)-old_bandt;

% useractive = zeros(I,T);
% for i = 1:I
%     useractive(i,timesplit(i,3):timesplit(i,4)) = 1;
% end
% 
% bandactive = zeros(J,T);
% for j = 1:J
%     bandactive(j,bandtime(j,2):bandtime(j,3)) = 1;
% end
            
% A(any(A==5, 2),:)=[]            % method to remove rows    

% test
% find(friend == 2 & band == 61)
% test2
% A =  find(friend == 2 & band == 555);
% b = find(friend == 2 & band == 555 & timeobs == 318)
% DV(b)
% DV(b+1)
% band(b+1)
% timeobs(b+1)

% member_a = [];
% friend_a = [];
% band_a = [];
% % timeobs = [];
% DV1 = [];
% largeNum_a = 48000000
largeNum_r = 900000

friendlist_updated = [];

% member_r = zeros(largeNum_r,1);
% friend_r = zeros(largeNum_r,1);
% band_r = zeros(largeNum_r,1);
% timeobs = zeros(largeNum_r,1);
% DV = zeros(largeNum_r,1);
% % obs_start = zeros(size(friendlist,1),J);
% friend_adopted = zeros(largeNum_r,1);
% % bandpopular_a = zeros(largeNum_r,1);
% totaladopt_r = zeros(largeNum_r,1);
% sinceintro_r = zeros(largeNum_r,1);
% explorerm_r = zeros(largeNum_r,1);
% explorerf_r = zeros(largeNum_r,1);
% EAm_r = zeros(largeNum_r,1);
% MAm_r = zeros(largeNum_r,1);
% LAm_r = zeros(largeNum_r,1);
% EAf_r = zeros(largeNum_r,1);
% MAf_r = zeros(largeNum_r,1);
% LAf_r = zeros(largeNum_r,1);

% ind = 1;
for i = 1:size(friendlist,1)                                         % i here  is the row num of friendlist
    member_id = friendlist(i,1);
    friend_id = friendlist(i,2);
    
    member_r = zeros(largeNum_r,1);
    friend_r = zeros(largeNum_r,1);
    band_r = zeros(largeNum_r,1);
    timeobs = zeros(largeNum_r,1);
    DV = zeros(largeNum_r,1);
    friend_adopted = zeros(largeNum_r,1);
    totaladopt_r = zeros(largeNum_r,1);
    sinceintro_r = zeros(largeNum_r,1);
    explorerm_r = zeros(largeNum_r,1);
    explorerf_r = zeros(largeNum_r,1);
    EAm_r = zeros(largeNum_r,1);
    MAm_r = zeros(largeNum_r,1);
    LAm_r = zeros(largeNum_r,1);
    EAf_r = zeros(largeNum_r,1);
    MAf_r = zeros(largeNum_r,1);
    LAf_r = zeros(largeNum_r,1);
    
    ind = 1;
    for j = 1:J
        %         row_adoptfij = adoptions(:,1) == id_friend & adoptions(:,2) == j;
        %         adoptfij_time = adoptions(row_adoptfij,3);
        adoptfij_time = band_adopt_mat(friend_id,j);
        pre_start = max([timesplit(member_id,3) timesplit(friend_id,3) bandtime(j,2)]);
        if adoptfij_time~=0
            pre_end = min([timesplit(member_id,4) adoptfij_time bandtime(j,3)]);
        else
            pre_end = min([timesplit(member_id,4) timesplit(friend_id,4) bandtime(j,3)]);
        end
        interval = pre_end - pre_start +1;
        if interval >= 1
%             if isempty(adoptfij_time)==0 && (adoptfij_time < pre_start || adoptfij_time > pre_end)      % remove obs of band j if it was adopted by i in period 1
            if adoptfij_time~=0 && (adoptfij_time < pre_start || adoptfij_time > pre_end)
            else
                %                 obs_start(i,j) = pre_start;
%                 member_a(ind) = id_member;
%                 friend_a(ind) = id_friend;
                member_r(ind:ind+interval-1) = member_id;
                friend_r(ind:ind+interval-1) = friend_id;
                band_r(ind:ind+interval-1) = j;
                timeobs(ind:ind+interval-1) = pre_start:pre_end;                
                if adoptfij_time~=0
%                     interval_r = adoptfij_time- pre_start +1;
%                 if isempty(adoptfij_time)==0                                                    % w: NEED adoptfij_time HERE !!!
%                     DV(ind+interval_r-1) = 1;
                    DV(ind+interval-1) = 1;
                    if ind+interval-1 > largeNum_r
                        DV(ind:ind+interval-2) = 0;
                    end
                else
                    if ind+interval-1 > largeNum_r
                        DV(ind:ind+interval-1) = 0;
                    end                        
                end                
                adoptmij_time = band_adopt_mat(member_id, j);                
                % if adoptmij_time~=0 && adoptmij_time < adoptfij_time && adoptmij_time > pre_start
                if adoptmij_time~=0 && adoptmij_time < pre_end && adoptmij_time > pre_start
                    %         rows_ind2 = member_t == member_id & friend_t == friend_id & band_t == band_id & timeobs >= adoptmij_time & timeobs - adoptmij_time <= 8;
                    rows_ind_r = timeobs(ind:ind+interval-1) >= adoptmij_time & timeobs(ind:ind+interval-1) - adoptmij_time <= 8;    % length(rows_ind_r) = length(ind:ind+interval-1)
                    if ind+interval-1 <= largeNum_r
                        rows_indall = logical([sparse(ind-1,1); rows_ind_r; sparse(largeNum_r-(ind+interval-1),1)]);
                    else
                        rows_indall = logical([sparse(ind-1,1); rows_ind_r; zeros(largeNum_r-(ind+interval-1),1)]);
                    end
                    friend_adopted(rows_indall) = 1;
                else
                    if ind+interval-1 > largeNum_r
                        friend_adopted(ind:ind+interval-1) = 0;
                    end
                end
                
                if pre_start > 1 && pre_end > 1
                    totaladopt_r(ind:ind+interval-1) = band_pop(pre_start-1:pre_end-1,j);                  % fixed bug here
                end
                sinceintro_r(ind:ind+interval-1) = timeobs(ind:ind+interval-1)-bandtime(j,2);
                
                explorerm_r(ind:ind+interval-1) = explorer(member_id);
                explorerf_r(ind:ind+interval-1) = explorer(friend_id);
                EAm_r(ind:ind+interval-1) = EAcat(member_id,1);
                MAm_r(ind:ind+interval-1) = EAcat(member_id,2);
                LAm_r(ind:ind+interval-1) = EAcat(member_id,3);
                EAf_r(ind:ind+interval-1) = EAcat(friend_id,1);
                MAf_r(ind:ind+interval-1) = EAcat(friend_id,2);
                LAf_r(ind:ind+interval-1) = EAcat(friend_id,3);
                
                ind = ind+interval;
                
            end
        else
        end
    end
    if ind > 1
        friendlist_updated = [friendlist_updated;[member_id, friend_id]];
    else
    end
    if ind <= largeNum_r
        member_r = member_r(1:(ind-1));
        friend_r = friend_r(1:(ind-1));
        band_r = band_r(1:(ind-1));
        timeobs = timeobs(1:(ind-1));
        DV = DV(1:(ind-1));
        friend_adopted = friend_adopted(1:(ind-1));

        totaladopt_r = totaladopt_r(1:(ind-1));
        sinceintro_r = sinceintro_r(1:(ind-1));
        explorerm_r = explorerm_r(1:(ind-1));
        explorerf_r = explorerf_r(1:(ind-1));
        EAm_r = EAm_r(1:(ind-1));
        MAm_r = MAm_r(1:(ind-1));
        LAm_r = LAm_r(1:(ind-1));
        EAf_r = EAf_r(1:(ind-1));
        MAf_r = MAf_r(1:(ind-1));
        LAf_r = LAf_r(1:(ind-1));
    else
    end
    
%     mat_r = [member_r friend_r band_r timeobs DV friend_adopted totaladopt_r sinceintro_r explorerm_r explorerf_r EAm_r MAm_r LAm_r EAf_r MAf_r LAf_r];
% %     save(sprintf('%d_%d',member_id, i), 'mat_r');
% %    save(sprintf('C:\\Users\\etp21998@eur.nl\\matfolder\\%d_%d.mat',member_id, i),'mat_r')      % save in a different folder
% %    save(sprintf('E:\\Wei\\Graduate\\Matlab\\matfolder\\%d_%d.mat',member_id, i),'mat_r') 
%     save(sprintf('E:\\matfolder\\%d_%d.mat',member_id, i),'mat_r')
    
    if ind > 1
        updated_i = size(friendlist_updated,1);
        mat_r = [member_r friend_r band_r timeobs DV friend_adopted totaladopt_r sinceintro_r explorerm_r explorerf_r EAm_r MAm_r LAm_r EAf_r MAf_r LAf_r];
        save(sprintf('E:\\matfolder\\%d_%d.mat',member_id,updated_i), 'mat_r');
        clearvars mat_r
    else
    end
    
    clearvars mat_r
    
    if i >=156
        break
    end
end

csvwrite('E:\matfolder\updated_friend.csv',friendlist_updated);
save('E:\matfolder\updated_friend','friendlist_updated')
% member_r = member_r(1:(ind-1));
% friend_r = friend_r(1:(ind-1));
% band_r = band_r(1:(ind-1));
% timeobs = timeobs(1:(ind-1));
% DV = DV(1:(ind-1));
% friend_adopted = friend_adopted(1:(ind-1));
% 
% totaladopt_r = totaladopt_r(1:(ind-1));
% sinceintro_r = sinceintro_r(1:(ind-1));
% explorerm_r = explorerm_r(1:(ind-1));
% explorerf_r = explorerf_r(1:(ind-1));
% EAm_r = EAm_r(1:(ind-1));
% MAm_r = MAm_r(1:(ind-1));
% LAm_r = LAm_r(1:(ind-1));
% EAf_r = EAf_r(1:(ind-1));
% MAf_r = MAf_r(1:(ind-1));
% LAf_r = LAf_r(1:(ind-1));

display('DV')


% test
find(sinceintro_r<0)


%%
% mat_r = [member_r friend_r band_r timeobs DV friend_adopted totaladopt_r sinceintro_r explorerm_r explorerf_r EAm_r MAm_r LAm_r EAf_r MAf_r LAf_r];
% % mat2 = [member_t friend_t band_t timeobs DV2 friend_adopted totaladopt_t sinceintro_t explorerm_t explorerf_t EAm_t MAm_t LAm_t EAf_t MAf_t LAf_t];
% 
% % clearvars -EXCEPT mat1 mat2
% 
% % save('mat1.mat','mat1')
% % save('mat2.mat','mat2')
% save('mat_r.mat','mat_r', '-v7.3') ;
% % save('mat2.mat','mat2', '-v7.3') ;

display('save as csv')
% clear all

% csvwrite('mat1.csv',mat1);
% csvwrite('mat2.csv',mat2);

%% try load saved files

% try save in different folder
% save('E:\matfolder\randomntrial','friendlist_updated')

% load('2898_1.mat');
% load('2898_2.mat');
% fri_row1 = load('2898_1.mat');               % fri_row1 is struct object
% fri_row2 = load('2898_2.mat');
% fri_row1.mat_r(1,2)
% fri_row2.mat_r(1,2)
% 
% 
% load(sprintf('%d_%d',2898,3))
% friend_files = cell(3,1);
% for i = 1:3
%     friend_files{i} = load(sprintf('C:\\Users\\etp21998@eur.nl\\matfolder\\%d_%d.mat',2898,i));
% end
% friend_files{1}.mat_r(1,2)


% whos(variablename)                      % see how much memory space used by the variable

% load('E:\Wei\Graduate\Matlab\matfolder\randomname.mat')

%% EA cols
% friendlist_updated = csvread('C:\Users\etp21998@eur.nl\matfolder\updated_friend.csv');
% member_list = csvread('mlist.csv');
% EAcat = csvread('EAdummy3.csv');
% 
% EA_location = cell(length(member_list),1);
% f_type = zeros(length(member_list),3);
% m_type = zeros(length(member_list),3);
% for i = 1:length(member_list)
%     member = member_list(i);
%     mf_rows = find(friendlist_updated(:,1)==member);
%     friends = friendlist_updated(mf_rows,2);
%     EAall = EAcat(friends,:);
%     EAsum = sum(EAall,1);
%     m_type(i,:) = EAcat(member,:);
%     f_EA = EAsum(1)~=0;
%     f_MA = EAsum(2)~=0;
%     f_LA = EAsum(3)~=0;
%     f_type(i,:) = [f_EA f_MA f_LA];
%     if f_EA == 1 && f_MA == 1 && f_LA == 1
%         EA_col = [14:15];
%     elseif f_EA == 1 && f_MA == 1
%         EA_col = 14;
%     elseif f_MA == 1 && f_LA == 1
%         EA_col = 15;
%     elseif f_EA == 1 && f_LA == 1
%         EA_col = 14;
%     else
%         EA_col = [];
%     end
%     EA_location{i} = EA_col; 
% end




