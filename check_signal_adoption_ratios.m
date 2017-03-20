load('matpstd2.mat');

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
for i= 2:size(matp,1)
    current_obs = matp(i,3);
    if current_obs~=last_obs || any(f_change_ind==i)
        f_adopt = [f_adopt sum(matp(start_ind:(i-1),5))];
        start_ind = i;
    else
    end
    last_obs = current_obs;
end
f_adopt = [f_adopt sum(matp(start_ind:i,5))];     

length(f_adopt)                   % number of adoption observations in DV of the data matrix used for modeling
sum(f_adopt==1)
sum(matp(:,5))

% check why obs only include those bands adopted by both friend and member

%% testing 
tic
A = 1:9
B = [1 3]
k = 0;
for i=1:length(B)
    k= k+any(A==B(i));
end
k
toc

tic
A = 1:9
B = [1 3]
C = repmat(B',1,length(A));
D = repmat(A,length(B),1);
sum(sum(C==D))
toc

%% test cell
mycell = cell(2,1)
mycell{1} = 1:9;
mycell{2} = [1 3];
mycell{2}(1,2)

tic
A = mycell{1};
B = mycell{2};
k = 0;
for i=1:length(B)
    k= k+any(A==B(i));
end
k
toc

%% non-strict rule of adoption (new bands not appeared in the first 104 weeks + adopted by member) using bandadoption3.csv data
% compute the ratio between witnessed new bands adopted by the member and adopted new bands by the friend
% adoption signal ratio

adoptions = csvread('bandadoptions3.csv');

last_obs = adoptions(1,1);
start_ind = 1;
stop_ind = [];
for i= 2:size(adoptions,1)
    current_obs = adoptions(i,1);
    if current_obs~=last_obs
        stop_ind = [stop_ind i-1]; 
        start_ind = [start_ind i];
    else
    end
    last_obs = current_obs;
end     
stop_ind = [stop_ind size(adoptions,1)];

I = 8320;
store_witnessed_bands = cell(I,1);
for i = 1:I
    store_witnessed_bands{i} = adoptions(start_ind(i):stop_ind(i),2);
end

friendlist = csvread('friends3.csv');
num_new_bands_adopt_by_member = zeros(size(friendlist,2),1);
num_new_bands_adopt_by_friend = zeros(size(friendlist,2),1);
num_new_bands_adopt_by_both_member_friend = zeros(size(friendlist,2),1);
for k = 1:size(friendlist,1)
    member_ind = friendlist(k,1);
    friend_ind = friendlist(k,2);
    num_new_bands_adopt_by_member(k) = length(store_witnessed_bands{member_ind});
    num_new_bands_adopt_by_friend(k) = length(store_witnessed_bands{friend_ind});
    overlap_new_bands = 0;
    for j = 1:length(store_witnessed_bands{friend_ind})
        overlap_new_bands = overlap_new_bands + any(store_witnessed_bands{member_ind}==store_witnessed_bands{friend_ind}(j));
    end
    num_new_bands_adopt_by_both_member_friend(k) = overlap_new_bands;
end
        
mean(num_new_bands_adopt_by_both_member_friend)
mean(num_new_bands_adopt_by_member)
mean(num_new_bands_adopt_by_friend)

witness_new_band_adoption_ratio = num_new_bands_adopt_by_both_member_friend./num_new_bands_adopt_by_member;
mean(witness_new_band_adoption_ratio)

%% strict rule of adoption (new bands appear after week 104 + adopted by member + not introduced before user listen period) using bandadoptions_strict_adopt.csv data

adoptions = csvread('bandadoptions_strict_adopt.csv',1,0); 

last_obs = adoptions(1,1);
start_ind = 1;
stop_ind = [];
for i= 2:size(adoptions,1)
    current_obs = adoptions(i,1);
    if current_obs~=last_obs
        stop_ind = [stop_ind i-1]; 
        start_ind = [start_ind i];
    else
    end
    last_obs = current_obs;
end     
stop_ind = [stop_ind size(adoptions,1)];

user_ids_7130 = adoptions(start_ind,1);

I = 8320;
store_witnessed_bands = cell(I,1);
for i = 1:I
    if isempty(find(user_ids_7130==i))
        store_witnessed_bands{i} = [];
    else
        store_position = find(user_ids_7130==i);
        store_witnessed_bands{i} = adoptions(start_ind(store_position):stop_ind(store_position),2);
    end
end

friendlist = csvread('new_friendlist_7623.csv');
num_new_bands_adopt_by_member = zeros(size(friendlist,2),1);
num_new_bands_adopt_by_friend = zeros(size(friendlist,2),1);
num_new_bands_adopt_by_both_member_friend = zeros(size(friendlist,2),1);
for k = 1:size(friendlist,1)
    member_ind = friendlist(k,1);
    friend_ind = friendlist(k,2);
    num_new_bands_adopt_by_member(k) = length(store_witnessed_bands{member_ind});
    num_new_bands_adopt_by_friend(k) = length(store_witnessed_bands{friend_ind});
    overlap_new_bands = 0;
    for j = 1:length(store_witnessed_bands{friend_ind})
        overlap_new_bands = overlap_new_bands + any(store_witnessed_bands{member_ind}==store_witnessed_bands{friend_ind}(j));
    end
    num_new_bands_adopt_by_both_member_friend(k) = overlap_new_bands;
end
        
mean(num_new_bands_adopt_by_both_member_friend)
mean(num_new_bands_adopt_by_member)
mean(num_new_bands_adopt_by_friend)

witness_new_band_adoption_ratio = num_new_bands_adopt_by_both_member_friend./num_new_bands_adopt_by_member;
mean(witness_new_band_adoption_ratio)

load('matp_trend_strict_rm.mat');
sum(matp(:,5))

%% less strict rule of adoption (new bands appear after week 104 + adopted by member + not introduced before user listen period + introduced before user listen period but user has not listened to the band in the 1st year of her listen history)
adoptions = csvread('bandadoptions_lenient_adopt.csv',1,0); 

last_obs = adoptions(1,1);
start_ind = 1;
stop_ind = [];
for i= 2:size(adoptions,1)
    current_obs = adoptions(i,1);
    if current_obs~=last_obs
        stop_ind = [stop_ind i-1]; 
        start_ind = [start_ind i];
    else
    end
    last_obs = current_obs;
end     
stop_ind = [stop_ind size(adoptions,1)];

user_ids_7481 = adoptions(start_ind,1);

I = 8320;
store_witnessed_bands = cell(I,1);
for i = 1:I
    if isempty(find(user_ids_7481==i))
        store_witnessed_bands{i} = [];
    else
        store_position = find(user_ids_7481==i);
        store_witnessed_bands{i} = adoptions(start_ind(store_position):stop_ind(store_position),2);
    end
end

friendlist = csvread('new_friendlist_8088.csv',1,0);
num_new_bands_adopt_by_member = zeros(size(friendlist,2),1);
num_new_bands_adopt_by_friend = zeros(size(friendlist,2),1);
num_new_bands_adopt_by_both_member_friend = zeros(size(friendlist,2),1);
for k = 1:size(friendlist,1)
    member_ind = friendlist(k,1);
    friend_ind = friendlist(k,2);
    num_new_bands_adopt_by_member(k) = length(store_witnessed_bands{member_ind});
    num_new_bands_adopt_by_friend(k) = length(store_witnessed_bands{friend_ind});
    overlap_new_bands = 0;
    for j = 1:length(store_witnessed_bands{friend_ind})
        overlap_new_bands = overlap_new_bands + any(store_witnessed_bands{member_ind}==store_witnessed_bands{friend_ind}(j));
    end
    num_new_bands_adopt_by_both_member_friend(k) = overlap_new_bands;
end
        
mean(num_new_bands_adopt_by_both_member_friend)
mean(num_new_bands_adopt_by_member)
mean(num_new_bands_adopt_by_friend)

witness_new_band_adoption_ratio = num_new_bands_adopt_by_both_member_friend./num_new_bands_adopt_by_member;
mean(witness_new_band_adoption_ratio)

