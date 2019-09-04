store_count = [];
store_ind = [];
ind = friendlist(1,1);
count = 0;
for i = 1:size(friendlist,1)
    prev_ind = ind;
    ind = friendlist(i,1);
    if prev_ind~=ind
        store_count = [store_count; count];
        store_ind = [store_ind; ind];
        count = 1;
    else
        count = count+1;
    end
end
store_count = [store_count; count];
store_ind = [store_ind; ind];

%%
load('matp_friend_reverse.mat');
predict_trend = csvread('predict_trend_log4061_lenient.csv',1,0);

index = find(ismember(matp(:,3),predict_trend(:,1)));
matp = matp(index,:);

load('indx.mat');
load('dummy_mat.mat');

matp = matp(indx,:);

%%
length(unique(matp(:,1)))
memebers = unique(matp(:,1));

temp= ismember(store_ind,memebers);
mb = store_ind(temp);
m_count= store_count(temp);

histogram(m_count,30)
xlabel('number of friends')
ylabel('count users')

%%
x1 = [0.00108	0.00123	0.00155]
x2 = [0.00090	0.00112	0.00163]
x3 = [0.00076	0.00102	0.00171]
y = [0.16	0.5	   0.84]
plot(y,x1)
hold on
plot(y,x2)
hold on
plot(y,x3)
hold off
legend('Influencee Innovativeness = 0.16','Influencee Innovativeness = 0.50','Influencee Innovativeness = 0.84','Location','northwest')
ylabel('band sampling probability') 
xlabel('Influencer Innovativeness percentile (0.16, 0.5, 0.84)') 
