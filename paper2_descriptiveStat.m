count_friends = csvread('count_friends8088_repeat.csv',1,0);

load('m_innov.mat');
load('f_innov.mat');

friendlist = csvread('new_friendlist_8088.csv',1,0);

scatter(m_innov,f_innov,'.')
xlabel('focal user innovativeness');
ylabel('friends innovativeness');

[R,P] = corrcoef(m_innov,f_innov)

combi = [friendlist(:,1) m_innov f_innov];

temp = combi(end,:);
combi(end,:) = [];
combi = [temp; combi];

m = zeros(size(count_friends,1),1);
mi = zeros(size(count_friends,1),1);
fi = zeros(size(count_friends,1),1);
r = 1;
ind = 1;
while r <= size(count_friends,1)
    m(r) = combi(ind,1);
    rows = count_friends(find(count_friends(:,1)==m(r)),2);
    mi(r) = combi(ind,2);
    fi(r) = mean(combi(ind:ind+rows-1,3));
    ind = ind+rows;
    r = r+1;
end
    
corr(count_friends(:,2),mi)
corr(count_friends(:,2),fi)

[R,P] = corrcoef(count_friends(:,2),mi)
[R,P] = corrcoef(count_friends(:,2),fi)

count_bands = csvread('bands_per_user_8320.csv',1,0);

mbands = zeros(size(friendlist,1),1);
fbands = zeros(size(friendlist,1),1);
for i = 1:size(friendlist,1)
    mbands(i) = count_bands(friendlist(i,1),2);
    fbands(i) = count_bands(friendlist(i,2),2);
end
corr(m_innov,mbands)
corr(f_innov,fbands)

[C,ia,ic] = unique(friendlist(:,1),'rows');
corr(m_innov(ia),mbands(ia))

[R,P] = corrcoef([m_innov(ia); f_innov],[mbands(ia); fbands])

