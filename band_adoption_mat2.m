clc
clear

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
largeNum_a = 600000

member_a = zeros(largeNum_a,1);
friend_a = zeros(largeNum_a,1);
band_a = zeros(largeNum_a,1);
% timeobs = [];
DV1 = zeros(largeNum_a,1);
% obs_start = zeros(size(friendlist,1),J);
% friend_adopted = [];
bandpopular_a = zeros(largeNum_a,1);
explorerm_a = zeros(largeNum_a,1);
explorerf_a = zeros(largeNum_a,1);
EAm_a = zeros(largeNum_a,1);
MAm_a = zeros(largeNum_a,1);
LAm_a = zeros(largeNum_a,1);
EAf_a = zeros(largeNum_a,1);
MAf_a = zeros(largeNum_a,1);
LAf_a = zeros(largeNum_a,1);

ind = 1;
for i = 1:size(friendlist,1)                                         % i here  is the row num of friendlist
    id_member = friendlist(i,1);
    id_friend = friendlist(i,2);
    for j = 1:J
        %         row_adoptfij = adoptions(:,1) == id_friend & adoptions(:,2) == j;
        %         adoptfij_time = adoptions(row_adoptfij,3);
        adoptfij_time = band_adopt_mat(id_friend,j);
        pre_start = max([timesplit(id_member,3) timesplit(id_friend,3) bandtime(j,2)]);
        pre_end = min([timesplit(id_member,4) timesplit(id_friend,4) bandtime(j,3)]);
        interval = pre_end - pre_start +1;
        if interval >= 1
%             if isempty(adoptfij_time)==0 && (adoptfij_time < pre_start || adoptfij_time > pre_end)      % remove obs of band j if it was adopted by i in period 1
            if adoptfij_time~=0 && (adoptfij_time < pre_start || adoptfij_time > pre_end)
            else
                %                 obs_start(i,j) = pre_start;
                member_a(ind) = id_member;
                friend_a(ind) = id_friend;
                band_a(ind) = j;
                if adoptfij_time~=0
%                 if isempty(adoptfij_time)==0                                                    % w: NEED adoptfij_time HERE !!!
                    DV1(ind) = 1;
                else
                end
                if pre_start > 1
                    bandpopular_a(ind) = cumu_diffusion(pre_start-1,j);
                else
                end
                explorerm_a(ind) = explorer(id_member);
                explorerf_a(ind) = explorer(id_friend);
                EAm_a(ind) = EAcat(id_member,1);
                MAm_a(ind) = EAcat(id_member,2);
                LAm_a(ind) = EAcat(id_member,3);
                EAf_a(ind) = EAcat(id_friend,1);
                MAf_a(ind) = EAcat(id_friend,2);
                LAf_a(ind) = EAcat(id_friend,3);
                ind = ind+1;
            end
        else
        end
    end
        if i >=55
            break
        end
end

member_a = member_a(1:(ind-1));
friend_a = friend_a(1:(ind-1));
band_a = band_a(1:(ind-1));
DV1 = DV1(1:(ind-1));
bandpopular_a = bandpopular_a(1:(ind-1));

explorerm_a = explorerm_a(1:(ind-1));
explorerf_a = explorerf_a(1:(ind-1));
EAm_a = EAm_a(1:(ind-1));
MAm_a = MAm_a(1:(ind-1));
LAm_a = LAm_a(1:(ind-1));
EAf_a = EAf_a(1:(ind-1));
MAf_a = MAf_a(1:(ind-1));
LAf_a = LAf_a(1:(ind-1));

display('DV1')
%%
% test: some member and friends have no overlap in obs period 2
% A = find(member == 66 & friend == 33);

% largeNum_t = 64000000
largeNum_t = 800000
member_t = zeros(largeNum_t,1);
friend_t = zeros(largeNum_t,1);
band_t = zeros(largeNum_t,1);
timeobs = zeros(largeNum_t,1);
DV2 = zeros(largeNum_t,1);
friend_adopted = zeros(largeNum_t,1);

totaladopt_t = zeros(largeNum_t,1);
sinceintro_t = zeros(largeNum_t,1);
explorerm_t = zeros(largeNum_t,1);
explorerf_t = zeros(largeNum_t,1);
EAm_t = zeros(largeNum_t,1);
MAm_t = zeros(largeNum_t,1);
LAm_t = zeros(largeNum_t,1);
EAf_t = zeros(largeNum_t,1);
MAf_t = zeros(largeNum_t,1);
LAf_t = zeros(largeNum_t,1);

adopted_bands = find(DV1);
ind = 1;
for r = 1:length(adopted_bands)
    row_ind = adopted_bands(r);
    member_id = member_a(row_ind);
    friend_id = friend_a(row_ind);
    band_id = band_a(row_ind);
    %     row_adoptfij = adoptions(:,1) == friend_id & adoptions(:,2) == band_id;
    %     adoptfij_time = adoptions(row_adoptfij,3);
    adoptfij_time = band_adopt_mat(friend_id, band_id);
    pre_start = max([timesplit(member_id,3) timesplit(friend_id,3) bandtime(band_id,2)]);
    pre_end = min([timesplit(member_id,4) adoptfij_time bandtime(band_id,3)]);
    interval = pre_end - pre_start+1;
    
    member_t(ind:ind+interval-1) = member_id;
    friend_t(ind:ind+interval-1) = friend_id;
    band_t(ind:ind+interval-1) = band_id;
    timeobs(ind:ind+interval-1) = pre_start:pre_end;
    DV2(ind+interval-1) = 1;                        % assign 1 to DV at adoption time t, get rid of obs after i adopts j
    
    %     row_adoptmij = adoptions(:,1) == member_id & adoptions(:,2) == band_id;
    %     adoptmij_time = adoptions(row_adoptmij,3);
    adoptmij_time = band_adopt_mat(member_id, band_id);
    %     if isempty(row_adoptmij)==0 && adoptfij_time > adoptmij_time      % w: ????
%      if isempty(adoptmij_time)==0 && adoptfij_time > adoptmij_time
    if adoptmij_time~=0 && adoptfij_time > adoptmij_time
%         rows_ind2 = member_t == member_id & friend_t == friend_id & band_t == band_id & timeobs >= adoptmij_time & timeobs - adoptmij_time <= 8;
        rows_ind2 = timeobs(ind:ind+interval-1) >= adoptmij_time & timeobs(ind:ind+interval-1) - adoptmij_time <= 8;
        rows_indall = logical([sparse(ind-1,1); rows_ind2; sparse(largeNum_t-(ind+interval-1),1)]);
        friend_adopted(rows_indall) = 1;
    else
    end
    if pre_start > 1 && pre_end > 1
        totaladopt_t(ind:ind+interval-1) = band_pop(pre_start-1:pre_end-1,band_id);                  % fixed bug here
    end
    sinceintro_t(ind:ind+interval-1) = timeobs(ind:ind+interval-1)-bandtime(band_id,2);
    
    explorerm_t(ind:ind+interval-1) = explorer(member_id);
    explorerf_t(ind:ind+interval-1) = explorer(friend_id);
    EAm_t(ind:ind+interval-1) = EAcat(member_id,1);
    MAm_t(ind:ind+interval-1) = EAcat(member_id,2);
    LAm_t(ind:ind+interval-1) = EAcat(member_id,3);
    EAf_t(ind:ind+interval-1) = EAcat(friend_id,1);
    MAf_t(ind:ind+interval-1) = EAcat(friend_id,2);
    LAf_t(ind:ind+interval-1) = EAcat(friend_id,3);
    
    ind = ind+interval;
end

member_t = member_t(1:(ind-1));
friend_t = friend_t(1:(ind-1));
band_t = band_t(1:(ind-1));
timeobs = timeobs(1:(ind-1));
DV2 = DV2(1:(ind-1));
friend_adopted = friend_adopted(1:(ind-1));

totaladopt_t = totaladopt_t(1:(ind-1));
sinceintro_t = sinceintro_t(1:(ind-1));
explorerm_t = explorerm_t(1:(ind-1));
explorerf_t = explorerf_t(1:(ind-1));
EAm_t = EAm_t(1:(ind-1));
MAm_t = MAm_t(1:(ind-1));
LAm_t = LAm_t(1:(ind-1));
EAf_t = EAf_t(1:(ind-1));
MAf_t = MAf_t(1:(ind-1));
LAf_t = LAf_t(1:(ind-1));

% test
% w = find(friend_adopted);
% friend_adopted(151)
% timeobs(151)
% band_t(151)
% member_t(151)
% friend_t(151)
% sum(DV2)

display('DV2')
%%
% explorerm_a = zeros(length(member_a),1);
% for i = 1:194                                % number of members here !!
%     m_ind = find(member_a == i);
%     explorerm_a(m_ind) = explorer(i);
% end
% 
% explorerf_a = zeros(length(friend_a),1);
% for i = 1:194                                % number of friends here !!
%     f_ind = find(friend_a == i);
%     explorerf_a(f_ind) = explorer(i);
% end
% 
% EAm_a = zeros(length(member_a),1);
% MAm_a = zeros(length(member_a),1);
% LAm_a = zeros(length(member_a),1);
% for i = 1:194                                % number of members here !!
%     m_ind = find(member_a == i);
%     EAm_a(m_ind) = EAcat(i,1);
%     MAm_a(m_ind) = EAcat(i,2);
%     LAm_a(m_ind) = EAcat(i,3);
% end
% 
% EAf_a = zeros(length(friend_a),1);
% MAf_a = zeros(length(friend_a),1);
% LAf_a = zeros(length(friend_a),1);
% for i = 1:194                                % number of friends here !!
%     f_ind = find(friend_a == i);
%     EAf_a(f_ind) = EAcat(i,1);
%     MAf_a(f_ind) = EAcat(i,2);
%     LAf_a(f_ind) = EAcat(i,3);
% end

% bandpopular_a = zeros(length(band_a),1);            % revise here !!
% for i = 1:size(friendlist,1)                % i here is the row num of friendlist
%     for j = 1:J
%         if obs_start(i,j) > 0
%             mij_ind = find(member_a == friendlist(i,1) & friend_a == friendlist(i,2) & band_a == j);
%             if obs_start(i,j) > 1
%                 bandpopular_a(mij_ind) = cumu_diffusion(obs_start(i,j)-1,j);
%             else
%             end
%         else
%         end
%     end
% end


% totaladopt_t = zeros(length(timeobs),1);
% for j = 1:J
%     for t = 1:T
%         jt_ind = find(band_t == j & timeobs == t);
%         if t > 1
%             totaladopt_t(jt_ind) = band_pop(t-1,j);
%         else
%         end
%     end
% end
% 
% sinceintro_t = zeros(length(timeobs),1);
% for j = 1:J
%     j_ind = find(band_t == j);
%     sinceintro_t(j_ind) = timeobs(j_ind)-bandtime(j,2);
% end

% test
find(sinceintro_t<0)

% explorerm_t = zeros(length(member_t),1);
% for i = 1:194                                % number of members here !!
%     m_ind = find(member_t == i);
%     explorerm_t(m_ind) = explorer(i);
% end
% 
% explorerf_t = zeros(length(friend_t),1);
% for i = 1:194                                % number of friends here !!
%     f_ind = find(friend_t == i);
%     explorerf_t(f_ind) = explorer(i);
% end
% 
% EAm_t = zeros(length(member_t),1);
% MAm_t = zeros(length(member_t),1);
% LAm_t = zeros(length(member_t),1);
% for i = 1:194                                % number of members here !!
%     m_ind = find(member_t == i);
%     EAm_t(m_ind) = EAcat(i,1);
%     MAm_t(m_ind) = EAcat(i,2);
%     LAm_t(m_ind) = EAcat(i,3);
% end
% 
% EAf_t = zeros(length(friend_t),1);
% MAf_t = zeros(length(friend_t),1);
% LAf_t = zeros(length(friend_t),1);
% for i = 1:194                                % number of friends here !!
%     f_ind = find(friend_t == i);
%     EAf_t(f_ind) = EAcat(i,1);
%     MAf_t(f_ind) = EAcat(i,2);
%     LAf_t(f_ind) = EAcat(i,3);
% end

%%
mat1 = [member_a friend_a band_a DV1 bandpopular_a explorerm_a explorerf_a EAm_a MAm_a LAm_a EAf_a MAf_a LAf_a];
mat2 = [member_t friend_t band_t timeobs DV2 friend_adopted totaladopt_t sinceintro_t explorerm_t explorerf_t EAm_t MAm_t LAm_t EAf_t MAf_t LAf_t];

% clearvars -EXCEPT mat1 mat2

% save('mat1.mat','mat1')
% save('mat2.mat','mat2')
save('mat1.mat','mat1', '-v7.3') ;
save('mat2.mat','mat2', '-v7.3') ;

display('save as csv')
csvwrite('mat1.csv',mat1);
csvwrite('mat2.csv',mat2);

