clc
clear

adoptions = csvread('bandadoptions3.csv');
adopt_ones = ones(size(adoptions, 1),1);

% total number of weeks: length(105:527) = 423
T = 423
I = 8320
J = 6046
old_bandt = 104
band_adoption = cell(T,1);

for t = 1:T
    [r c v] = find(adoptions(:,3)==(t+old_bandt));    
    adopt_ind = adoptions(r,:);
    band_adoption{t} = sparse(adopt_ind(:,1),adopt_ind(:,2),adopt_ones(r),I,J);
end

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

explorer = zeros(I,1);
for i = 1:I
    for t = old_bandt+1:timesplit(i,3)             
        explorer(i) = explorer(i)+sum(band_adoption{t-old_bandt}(i,:));
    end
end

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

moving_avg = 32;
moving_w = zeros(J,T-moving_avg);
for j = 1:J
    for t = 1:T-moving_avg
        moving_w(j,t) = mean(diffusion_jt(t:t+moving_avg-1,j));
    end
end
sum(diffusion_jt(:,1577))                % Ed Sheeran
plot(moving_w(1577,:))
sum(diffusion_jt(:,799))                % Bruno Mars
plot(moving_w(799,:))

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

adoptingweek = csvread('adoptingweek3.csv');

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

useractive = zeros(I,T);
for i = 1:I
    useractive(i,timesplit(i,3):timesplit(i,4)) = 1;
end

bandactive = zeros(J,T);
for j = 1:J
    bandactive(j,bandtime(j,2):bandtime(j,3)) = 1;
end
            
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


member_a = [];
friend_a = [];
band_a = [];
% timeobs = [];
DV1 = [];
% friend_adopted = [];
for i = 1:size(friendlist,1)                                         % i here  is the row num of friendlist
    id_member = friendlist(i,1);
    id_friend = friendlist(i,2);
    for j = 1:J
        row_adoptfij = find(adoptions(:,1) == id_friend & adoptions(:,2) == j);    
        adoptfij_time = adoptions(row_adoptfij,3);
        if isempty(row_adoptfij)==0 && adoptfij_time < timesplit(id_friend,3)      % remove obs of band j if it was adopted by i in period 1
        else
            product = useractive(id_member,:).*useractive(id_friend,:).*bandactive(j,:);          % only include bands with the 3 overlaps
            if sum(product)>0
                member_a = [member_a; id_member];
                friend_a = [friend_a; id_friend];
                band_a = [band_a; j];
                %                 timeobs = [timeobs; find(product)'];
                DV1 = [DV1; 0];
                %                 friend_adopted = [friend_adopted; zeros(sum(product),1)];
                if isempty(row_adoptfij)==0                                            
                    row_ij = find(member_a == id_member & friend_a == id_friend & band_a == j);   % assign 1 to DV for band j
                    DV1(row_ij) = 1;
                    %                 rows_ind = find(member == id_member & friend == id_friend & band == j & timeobs > adoptfij_time);   % get rid of obs after i adopts j
                    %                 DV(rows_ind) = [];
                    %                 friend(rows_ind) = [];
                    %                 band(rows_ind) = [];
                    %                 timeobs(rows_ind) = [];
                    %                 member(rows_ind) = [];
                    %                 friend_adopted(rows_ind) = [];
                else
                end
                %             row_adoptmij = find(adoptions(:,1) == id_member & adoptions(:,2) == j);
                %             adoptmij_time = adoptions(row_adoptmij,3);
                %             if isempty(row_adoptmij)==0
                %                 rows_ind2 = find(member == id_member & friend == id_friend & band == j & timeobs - adoptmij_time <= 8);  % member == id_member ??
                %                 % rows_ind2 = find(friend == id_friend & band == j & timeobs - adoptmij_time <= 8);  % member == id_member ??
                %                 friend_adopted(rows_ind2) = 1;
                %             else
                %             end
            else
            end
        end
    end
    if i >=5
        break
    end
end

%%
% test: some member and friends have no overlap in obs period 2
% A = find(member == 66 & friend == 33);

member_t = [];
friend_t = [];
band_t = [];
timeobs = [];
DV2 = [];
friend_adopted = [];

adopted_bands = find(DV1);
for r = 1:length(adopted_bands)
    row_ind = adopted_bands(r);
    member_id = member_a(row_ind);
    friend_id = friend_a(row_ind);
    band_id = band_a(row_ind);
    product = useractive(member_id,:).*useractive(friend_id,:).*bandactive(band_id,:);
    member_t = [member_t; member_id.*ones(sum(product),1)];
    friend_t = [friend_t; friend_id.*ones(sum(product),1)];
    band_t = [band_t; band_id.*ones(sum(product),1)];
    timeobs = [timeobs; find(product)'];
    DV2 = [DV2; zeros(sum(product),1)];
    friend_adopted = [friend_adopted; zeros(sum(product),1)];
    
    row_adoptfij = find(adoptions(:,1) == friend_id & adoptions(:,2) == band_id);
    adoptfij_time = adoptions(row_adoptfij,3);
    row_ijt = find(member_t == member_id & friend_t == friend_id & band_t == band_id & timeobs == adoptfij_time);   % assign 1 to DV at adoption time t
    DV2(row_ijt) = 1;
    rows_ind = find(member_t == member_id & friend_t == friend_id & band_t == band_id & timeobs > adoptfij_time);   % get rid of obs after i adopts j
    DV2(rows_ind) = [];
    friend_t(rows_ind) = [];
    band_t(rows_ind) = [];
    timeobs(rows_ind) = [];
    member_t(rows_ind) = [];
    friend_adopted(rows_ind) = [];
    
    row_adoptmij = find(adoptions(:,1) == member_id & adoptions(:,2) == band_id);
    adoptmij_time = adoptions(row_adoptmij,3);
    if isempty(row_adoptmij)==0 && adoptfij_time > adoptmij_time           % w: adoptfij_time >= adoptmij_time ??
        rows_ind2 = find(member_t == member_id & friend_t == friend_id & band_t == band_id & timeobs >= adoptmij_time & timeobs - adoptmij_time <= 8);
        % rows_ind2 = find(friend == id_friend & band == j & timeobs - adoptmij_time <= 8);  % member == id_member ??
        friend_adopted(rows_ind2) = 1;
    else
    end
end

% test
% w = find(friend_adopted);
% friend_adopted(151)
% timeobs(151)
% band_t(151)
% member_t(151)
% friend_t(151)
% sum(DV2)

%%
explorerm_a = zeros(length(member_a),1);
for i = 1:194                                % number of members here !!
    m_ind = find(member_a == i);
    explorerm_a(m_ind) = explorer(i);
end

explorerf_a = zeros(length(friend_a),1);
for i = 1:194                                % number of friends here !!
    f_ind = find(friend_a == i);
    explorerf_a(f_ind) = explorer(i);
end

EAm_a = zeros(length(member_a),1);
MAm_a = zeros(length(member_a),1);
LAm_a = zeros(length(member_a),1);
for i = 1:194                                % number of members here !!
    m_ind = find(member_a == i);
    EAm_a(m_ind) = EAcat(i,1);
    MAm_a(m_ind) = EAcat(i,2);
    LAm_a(m_ind) = EAcat(i,3);
end

EAf_a = zeros(length(friend_a),1);
MAf_a = zeros(length(friend_a),1);
LAf_a = zeros(length(friend_a),1);
for i = 1:194                                % number of friends here !!
    f_ind = find(friend_a == i);
    EAf_a(f_ind) = EAcat(i,1);
    MAf_a(f_ind) = EAcat(i,2);
    LAf_a(f_ind) = EAcat(i,3);
end

totaladopt_t = zeros(length(timeobs),1);
for j = 1:J
    for t = 1:T
        jt_ind = find(band_t == j & timeobs == t);
        totaladopt_t(jt_ind) = band_pop(t,j);
    end
end

sinceintro_t = zeros(length(timeobs),1);
for j = 1:J
    j_ind = find(band_t == j);
    sinceintro_t(j_ind) = timeobs(j_ind)-bandtime(j,2);
end

% test
find(sinceintro_t<0)

explorerm_t = zeros(length(member_t),1);
for i = 1:194                                % number of members here !!
    m_ind = find(member_t == i);
    explorerm_t(m_ind) = explorer(i);
end

explorerf_t = zeros(length(friend_t),1);
for i = 1:194                                % number of friends here !!
    f_ind = find(friend_t == i);
    explorerf_t(f_ind) = explorer(i);
end

EAm_t = zeros(length(member_t),1);
MAm_t = zeros(length(member_t),1);
LAm_t = zeros(length(member_t),1);
for i = 1:194                                % number of members here !!
    m_ind = find(member_t == i);
    EAm_t(m_ind) = EAcat(i,1);
    MAm_t(m_ind) = EAcat(i,2);
    LAm_t(m_ind) = EAcat(i,3);
end

EAf_t = zeros(length(friend_t),1);
MAf_t = zeros(length(friend_t),1);
LAf_t = zeros(length(friend_t),1);
for i = 1:194                                % number of friends here !!
    f_ind = find(friend_t == i);
    EAf_t(f_ind) = EAcat(i,1);
    MAf_t(f_ind) = EAcat(i,2);
    LAf_t(f_ind) = EAcat(i,3);
end

