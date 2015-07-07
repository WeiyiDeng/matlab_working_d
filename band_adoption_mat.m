clc
clear

adoptions = csvread('bandadoptions.csv');
adopt_ones = ones(size(adoptions, 1),1);

% total number of weeks: length(105:527) = 423
T = 423
I = 194
J = 7600
old_bandt = 104
band_adoption = cell(T,1);

for t = 1:T
    [r c v] = find(adoptions(:,3)==(t+old_bandt));    
    adopt_ind = adoptions(r,:);
    band_adoption{t} = sparse(adopt_ind(:,1),adopt_ind(:,2),adopt_ones(r),I,J);
end

% test
sum_a = 0;
for t = 1:T
    sum_a = sum_a + sum(sum([band_adoption{t}]));
end
sum_a

% cwd = pwd;
% cd(tempdir);
% pack                % pack can only be used in the command line ??
% cd(cwd)

for t = 1:T
    diffusion_jt(t,:) = sum(band_adoption{t},1);           % T*J
end

% [ri,cj,vs] = find(diffusion_jt);
% [sm,sn] = size(diffusion_jt);
% diffusion_jt = sparse(ri,cj,vs,sm,sn);

cumu_diffusion = cumsum(diffusion_jt);                    % T*J

% test: 
% max(max(cumu_diffusion))
% [row,col]=find(cumu_diffusion==108)     % MGMT BAND_ID=4096
% cumu_diffusion(:,150)                   % Adele
Timew = 1:T;
scatter(Timew,cumu_diffusion(:,4096));
scatter(Timew,cumu_diffusion(:,150));
%%
timesplit = csvread('tsplit.csv');

% test: these people should already be those with new band listens
find(timesplit(:,2:4)<=old_bandt)

explorer = zeros(I,1);
for i = 1:I
    for t = old_bandt+1:timesplit(i,3)             
        explorer(i) = explorer(i)+sum(band_adoption{t-old_bandt}(i,:));
    end
end

% test: moving average of 8 weeks
Timeww = 1:T-8;
for t = 1:T-8
    wtry(t) = mean(diffusion_jt(t:t+7,150));
end
plot(wtry)
scatter(Timeww,wtry)
pk = find(wtry==max(wtry))
pkw = pk(1)+3

% 
moving_avg = 8;
peakj = zeros(J,1);
for j = 1:J
    moving_w = zeros(1,T-moving_avg);
    for t = 1:T-moving_avg
        moving_w(t) = mean(diffusion_jt(t:t+moving_avg-1,j));
    end
    peakw = find(moving_w==max(moving_w));
    if length(peakw) ==1
    peakj(j) = peakw + moving_avg/2-1;     % there are still multiple max peak adoption weeks, how to fix?
    else
        ind_pkw = ceil(length(peakw)/2);
        peakj(j) = peakw(ind_pkw) + moving_avg/2;         %-1;
    end
end

introdate = csvread('introdate.csv');

introdate(:,2) = introdate(:,2)-old_bandt;    

% test: some have same peak week and introweek
% E: J = 139, use diffusion_jt(:,139) to see week adoptions
% introdate = peakw = 241
testA = peakj-introdate(:,2);
find(testA<0)

peakj = [introdate(:,1) peakj];
adoptions(:,3) = adoptions(:,3)-old_bandt;
csvwrite('peak.csv',peakj);
csvwrite('introdate.csv',introdate);
csvwrite('mod_adoptions.csv',adoptions);
% give headers "BAND_ID", "USER_ID", "week_mod" etc.
% to the csv files in Emeditor before importing to SQL 

adoptingweek = csvread('adoptingweek.csv');

EAij = 1-(adoptingweek(:,3)-adoptingweek(:,4))./(adoptingweek(:,5)-adoptingweek(:,4));
% some have the same intro week and peak week, get NAN for 0 in denominator
% test: 
testC = adoptingweek(:,5)==adoptingweek(:,4);
find(testC)

ind_beforepk = EAij > 0;                  % also get rid of NANs here
ind_beforepk = logical(1-ind_beforepk);
EAij(ind_beforepk) = 0;

EAi = zeros(I,1);
for i = 1:I
    ind_i = find(adoptingweek(:,1)==i);
    EAi(i) = median(EAij(ind_i));
end

EAi = [(1:I)' EAi];
csvwrite('EAi.csv',EAi);