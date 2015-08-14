%%
% test: some member and friends have no overlap in obs period 2
% A = find(member == 66 & friend == 33);

largeNum_t = 64000000
% largeNum_t = 800000
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
%     if adoptmij_time~=0 && adoptfij_time > adoptmij_time
    if adoptmij_time~=0 && adoptmij_time < adoptfij_time && adoptmij_time > pre_start
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