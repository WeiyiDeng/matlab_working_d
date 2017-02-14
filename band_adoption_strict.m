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

% largeNum_p = 68000000
largeNum_p = 600000

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
        adoptmij_time = band_adopt_mat(member_id, j);
        if adoptfij_time~=0 && adoptmij_time~=0
            if adoptfij_time > timesplit(friend_id,3) && adoptmij_time > timesplit(member_id,3)
                pre_start = max([timesplit(friend_id,3) bandtime(j,2)]);
                pre_end = min([timesplit(friend_id,4) bandtime(j,3)]);        % still overlap period between friend obs and band obs ??
                interval = pre_end - pre_start +1;
                member_p(ind:ind+interval-1) = member_id;
                friend_p(ind:ind+interval-1) = friend_id;
                band_p(ind:ind+interval-1) = j;
                timeobs(ind:ind+interval-1) = pre_start:pre_end;
                DV(ind+(adoptfij_time-pre_start)) = 1;
                prob_adopt_week(ind:ind+interval-1) = prob_adoption_yearj(pre_start:pre_end,j);
                week_diff = (pre_start:pre_end) - adoptmij_time;
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
% abs_week_diff = abs_week_diff(1:(ind-1));
new_week_diff = new_week_diff(1:(ind-1));
A_week_ijt = A_week_ijt(1:(ind-1));

display('save as mat')

% matp = [member_p friend_p band_p timeobs DV prob_adopt_week abs_week_diff];
matp = [member_p friend_p band_p timeobs DV prob_adopt_week new_week_diff A_week_ijt];
% clearvars -EXCEPT matp
save('matp2.mat','matp', '-v7.3') ;
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
load('matp.mat');
row_ind = 1
row_fid = [];
row_interval = [];
new_row = matp(1,2);
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

%% week dummy
load('matp.mat')
week1_dummy = matp(:,7);
% sum(week1_dummy==1)
week1_dummy(week1_dummy~=1)=0;
matp(:,7) = week1_dummy;

%% similar distributions
hold off
y = [0.4322    0.3299    0.2496    0.1784    0.1469    0.1049    0.0893    0.0792]
x = [1:8].*5
plot(x,y)
hold on
x = 0:1:40;
plot(x,gampdf(x,0.9354,27.4738))
% plot(x,gampdf(x,1.1,30));
% plot(x,exppdf(x,30));
% plot(x,exppdf(x,10));
% plot(x,exppdf(x,50));
plot(x,exppdf(x,27.5));          % same as gampdf(x,1,27.5)
plot(x,gampdf(x,1,27.5));

b = 100
triang_distr = @(x) (b-x)*2/((b-1)*(b-1));
y = triang_distr(x);
y(x<1 | x>b)=0;
plot(x,y);

hold off
x = 1:10
y = [0.1742    0.1891    0.1542    0.1356    0.1285    0.1326    0.1384    0.1302    0.1151    0.1156]
plot(x,y)

% try second half data
x = 0:1:40;
plot(x,10.15*gampdf(x,0.9354,27.4738))
hold on
plot(x,5*gampdf(x,1.1928,exp(2.0272)))
hold off

%% continuous innovativeness and explorer score
load('matp2.mat');

innov = csvread('EAi3.csv');
explor = csvread('explorer3.csv');

innov = innov(:,2);

members_in_mat = unique(matp(:,1));
friends_in_mat = unique(matp(:,2));

innov_m = zeros(size(matp,1),1);
explor_m = zeros(size(matp,1),1);
for i = 1:length(members_in_mat)
    member_id = members_in_mat(i);
    innov_m(matp(:,1)==member_id) = innov(member_id);
    explor_m(matp(:,1)==member_id) = explor(member_id);
end

display('members EA continous')

innov_f = zeros(size(matp,1),1);
explor_f = zeros(size(matp,1),1);
for i = 1:length(friends_in_mat)
    friend_id = friends_in_mat(i);
    innov_f(matp(:,2)==friend_id) = innov(friend_id);
    explor_f(matp(:,2)==friend_id) = explor(friend_id);
end

display('friends EA continous')

clearvars matp

EA_continous = [innov_m innov_f explor_m explor_f];
save('EA_continous2.mat','EA_continous','-v7.3');
clearvars EA_continous2

innov_contin = [innov_m innov_m.^2 innov_f innov_f.^2 innov_m.*innov_f (innov_m.*innov_f).^2];
save('innov_contin2.mat','innov_contin','-v7.3');
clearvars innov_contin2

explor_contin = [explor_m explor_m.^2 explor_f explor_f.^2 explor_m.*explor_f (explor_m.*explor_f).^2];
save('explor_contin2.mat','explor_contin','-v7.3');

% %% standardize
% load('innov_contin.mat');
% innov_m = innov_contin(:,1);
% innov_f = innov_contin(:,3);
% innov_m = (innov_m-mean(innov_m))./std(innov_m);
% innov_f = (innov_f-mean(innov_f))./std(innov_f);
% innov_contin = [innov_m innov_m.^2 innov_f innov_f.^2 innov_m.*innov_f (innov_m.*innov_f).^2];
% save('innov_contin.mat','innov_contin','-v7.3');
% clearvars innov_contin
% 
% load('explor_contin.mat');
% explor_m = explor_contin(:,1);
% explor_f = explor_contin(:,3);
% explor_m = (explor_m-mean(explor_m))./std(explor_m);
% explor_f = (explor_f-mean(explor_f))./std(explor_f);
% explor_contin = [explor_m explor_m.^2 explor_f explor_f.^2 explor_m.*explor_f (explor_m.*explor_f).^2];
% save('explor_contin.mat','explor_contin','-v7.3');

%% cut into chunks
clear all
load('matp.mat');
matp = matp(:,5:7);
I = size(matp,1);
numc = 4;                        % number of chunks-1
lenc = fix(I/numc);       % length of each chunk
% C = zeros(lenc,size(innov_contin,2));
for i = 1:numc-1
    C_matp = matp((i-1)*lenc+1:i*lenc,:);
    save(sprintf('E:\\Wei\\Graduate\\Matlab\\matfolder\\matp_%d.mat',i), 'C_matp');
    clearvars C_matp
end
if rem(I,numc)>0
    C_matp = matp(i*lenc+1:end,:);
    save(sprintf('E:\\Wei\\Graduate\\Matlab\\matfolder\\matp_%d.mat',numc), 'C_matp');
else
end
clearvars C_matp matp

load('innov_contin.mat');
% I = size(innov_contin,1);
% numc = 3;                        % number of chunks-1
% lenc = fix(I/numc);       % length of each chunk
% C = zeros(lenc,size(innov_contin,2));
for i = 1:numc-1
    C_innov = innov_contin((i-1)*lenc+1:i*lenc,:);
    save(sprintf('E:\\Wei\\Graduate\\Matlab\\matfolder\\innov_contin_%d.mat',i), 'C_innov');
    clearvars C_innov
end
if rem(I,numc)>0
    C_innov = innov_contin(i*lenc+1:end,:);
    save(sprintf('E:\\Wei\\Graduate\\Matlab\\matfolder\\innov_contin_%d.mat',numc), 'C_innov');
else
end
clearvars C_innov innov_contin

load('explor_contin.mat');
% I = size(explor_contin,1);
% numc = 3;                        % number of chunks-1
% lenc = fix(I/numc);       % length of each chunk
% C = zeros(lenc,size(innov_contin,2));
for i = 1:numc-1
    C_explor = explor_contin((i-1)*lenc+1:i*lenc,:);
    save(sprintf('E:\\Wei\\Graduate\\Matlab\\matfolder\\explor_contin_%d.mat',i), 'C_explor');
    clearvars C_explor
end
if rem(I,numc)>0
    C_explor = explor_contin(i*lenc+1:end,:);
    save(sprintf('E:\\Wei\\Graduate\\Matlab\\matfolder\\explor_contin_%d.mat',numc), 'C_explor');
else
end
clearvars C_explor explor_contin

%% VIF
% load('matp.mat');
load('innov_contin.mat');
load('explor_contin.mat');
% X = matp(:,6:7);
X = [explor_contin(:,1:4) innov_contin(:,1:4)];
R = corrcoef(X);
VIF = diag(inv(R))'

innov = csvread('EAi3.csv');
explor = csvread('explorer3.csv');
scatter(innov(:,2),explor)
corr(innov(:,2),explor)

load('members_in_mat.mat');

innov_cut = innov(:,2)
sth = members_in_mat;
mytry1 = explor(sth);
mytry2 = innov_cut(sth);
corr(mytry1,mytry2)

A = rand(100,1).*100;
B = rand(100,1);
corr(A,B)
vecA = [];
vecB = [];
for i = 1:length(A)
    rtimes = round(rand(1)*10000);
    vecA = [vecA; repmat(A(i),rtimes,1)];
    vecB = [vecB; repmat(B(i),rtimes,1)];
end
corr(vecA,vecB)

%% standardize
load('matp2.mat');
matp(:,6) = (matp(:,6)-mean(matp(:,6)))./std(matp(:,6));
save('matpstd2.mat','matp','-v7.3');
clearvars matp

load('innov_contin2.mat');
innov_m = innov_contin(:,1);
innov_f = innov_contin(:,3);
% innov_m = (innov_m-mean(innov_m))./std(innov_m);
% innov_f = (innov_f-mean(innov_f))./std(innov_f);
innov_contin = [innov_m innov_m.^2 innov_f innov_f.^2 innov_m.*innov_f (innov_m.*innov_f).^2];
for i = 1:size(innov_contin,2)
    innov_contin(:,i) = (innov_contin(:,i)-mean(innov_contin(:,i)))./std(innov_contin(:,i));
end
save('innov_contin_std2.mat','innov_contin','-v7.3');
clearvars innov_contin

load('explor_contin2.mat');
explor_m = explor_contin(:,1);
explor_f = explor_contin(:,3);
% explor_m = (explor_m-mean(explor_m))./std(explor_m);
% explor_f = (explor_f-mean(explor_f))./std(explor_f);
explor_contin = [explor_m explor_m.^2 explor_f explor_f.^2 explor_m.*explor_f (explor_m.*explor_f).^2];
for i = 1:size(explor_contin,2)
    explor_contin(:,i) = (explor_contin(:,i)-mean(explor_contin(:,i)))./std(explor_contin(:,i));
end
save('explor_contin_std2.mat','explor_contin','-v7.3');

%% plot estimated main effect and quadratic terms coefficients
load('innov_contin_std2.mat');
mean(innov_contin)
std(innov_contin)
max(innov_contin)
min(innov_contin)
load('explor_contin_std2.mat');
mean(explor_contin)
std(explor_contin)
max(explor_contin)
min(explor_contin)

%
x = -1.77:0.01:3.58;
plot(x,-0.0105*x+0.0666*x.^2,'r')
hold on
x = -1.73:0.01:10.0;
plot(x, 0.1276*x+0.0555*x.^2,'y')
x = -0.92:0.01:1.9;
plot(x,0*x-0.0653*x.^2,'b')
x = -0.94:0.01:7.52;
plot(x,-0.0226*x+0.0052*x.^2,'g')
% plot(x,-0.0357*x+0.0066*x.^2,'g')      
hold off

%% try load partial varibles from matfile
load('innov_contin_std.mat');
p1 = innov_contin(:,1:2);
p2 = innov_contin(:,3:4);
p3 = innov_contin(:,5:6);
save innov_contin_std_3ps.mat p1 p2 p3 -v7.3;
clear innov_contin p1 p2 p3

load('explor_contin_std.mat');
p1 = explor_contin(:,1:2);
p2 = explor_contin(:,3:4);
p3 = explor_contin(:,5:6);
save explor_contin_std_3ps.mat p1 p2 p3 -v7.3;
clear explor_contin p1 p2 p3

vars = whos('-file','innov_contin_std_3ps.mat');
sth = load('innov_contin_std_3ps.mat', vars(2).name);
SNames = fieldnames(sth);
sth = getfield(sth, SNames{1});

load('innov_contin_std_3ps.mat', vars(2).name)
innov_contin = matfile('innov_contin_std_3ps.mat');
innov_ms = innov_contin.p1(:,1:2);

%% memmapfile
load('matp.mat');

fileID = fopen('mytry.dat','w');
fwrite(fileID, matp,'double');
fclose(fileID);

% m = memmapfile('mytry.dat','Format', 'double');   % w: will load the whole data matrix as a vector
% length(m.Data(:,1))

% m = memmapfile('mytry.dat', 'Format', {'double', [308581 8], 'im'});
% mydt = m.Data.im;
% length(m.Data.im(:,1))

mm = memmapfile('mytry.dat', 'Format', {'double', [308581 1], 'mj'},'Repeat', 8);
A = mm.Data(1).mj;          % first col
dtmean = [];
for i = 1:8
    dtmean = [dtmean mean(mm.Data(i).mj)];
end

%%
load('matp2.mat');
x = 0:max(matp(:,7));
y = zeros(length(x),1);
y(1)= sum(matp(:,7)==0);
for i = 1:max(matp(:,7))
    y(i+1) = sum(matp(:,7)==i);
end
bar(x,y)
