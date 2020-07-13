% robustness check
load('m_innov.mat');
load('f_innov.mat');
friendlist8088 = csvread('new_friendlist_8088.csv',1,0);

temp_m = m_innov(end);
m_innov(end) = [];
m_innov = [temp_m; m_innov];

temp_f = f_innov(end);
f_innov(end) = [];
f_innov = [temp_f; f_innov];

m_list = m_innov(1);
m_prev = m_innov(1);
m_start = 1;
m_end = [];
for i = 2:length(m_innov)
    m_present = m_innov(i);
    if m_prev ~= m_present
        m_list = [m_list; m_present];
        m_start = [m_start; i];
        m_end = [m_end; i-1];
    end
    m_prev = m_present;
end
m_end = [m_end; i];

f_innov_mean = zeros(length(m_list),1);
for m = 1:length(m_list)
    f_innov_mean(m) = mean(f_innov(m_start(m):m_end(m)));
end

m = histc(m_list,0:0.05:1);
f = histc(f_innov_mean,0:0.05:1);
plot(0:0.05:1,m)
hold on
plot(0:0.05:1,f)
hold on
hist(m_list)
hold on
hist(f_innov_mean)
hold off

% plot mean friend innov vs member innov
x = m_list;
y = f_innov_mean;
binRange = 0:0.05:1;
hcx = histcounts(x,[binRange Inf]);
hcy = histcounts(y,[binRange Inf]);
figure
bar(binRange,[hcx;hcy]')
xlabel('innovativeness scores') 
ylabel('counts') 
legend({'innov of members','mean innov of friends per member'},'Location','northeast')
