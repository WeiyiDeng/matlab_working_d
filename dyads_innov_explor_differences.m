load('matp2.mat');
% matp = matp(1:54921,:);           % try with small sample

innov = csvread('EAi3.csv');                % unstandardized
explor = csvread('explorer3.csv');          % unstandardized

innov = innov(:,2);

m_col = zeros(10000,1);
f_col = zeros(10000,1);
m_col(1) = matp(1,1);
f_col(1) = matp(1,2);
row_ind = 1;
for i = 2:size(matp,1)
    if matp(i,2) == f_col(row_ind)
    else
        row_ind = row_ind+1;
        m_col(row_ind) = matp(i,1);
        f_col(row_ind) = matp(i,2);
    end
end

m_col = m_col(1:row_ind,:);
f_col = f_col(1:row_ind,:);
mf_dyads = [m_col f_col];
save('mf_dyads.mat','mf_dyads','-v7.3');

m_innov = zeros(size(mf_dyads,1),1);
f_innov = zeros(size(mf_dyads,1),1);
m_explor = zeros(size(mf_dyads,1),1);
f_explor = zeros(size(mf_dyads,1),1);

for j = 1:size(mf_dyads,1)
    m_innov(j) = innov(mf_dyads(j,1));
    f_innov(j) = innov(mf_dyads(j,2));
    m_explor(j) = explor(mf_dyads(j,1));
    f_explor(j) = explor(mf_dyads(j,2));
end

abs_diff_innov = abs(m_innov-f_innov);
abs_diff_explor = abs(m_explor-f_explor);

m_list = unique(mf_dyads(:,1));
innov_mode = zeros(length(m_list),1);
explor_mode = zeros(length(m_list),1);
innov_std = zeros(length(m_list),1);
explor_std = zeros(length(m_list),1);
innov_mean = zeros(length(m_list),1);
explor_mean = zeros(length(m_list),1);
for q = 1:length(m_list)
    innov_mode(q) = mode(abs_diff_innov(mf_dyads(:,1)==m_list(q)));
    innov_std(q) = std(abs_diff_innov(mf_dyads(:,1)==m_list(q)));
    innov_mean(q) = mean(abs_diff_innov(mf_dyads(:,1)==m_list(q)));
    explor_mode(q) = mode(abs_diff_explor(mf_dyads(:,1)==m_list(q)));
    explor_std(q) = std(abs_diff_explor(mf_dyads(:,1)==m_list(q)));
    explor_mean(q) = mean(abs_diff_explor(mf_dyads(:,1)==m_list(q)));
end

hist(innov_mean,50);
ylabel('num of dyads')
xlabel('mean of ? innovativeness')
hist(innov_std,50);
ylabel('num of dyads')
xlabel('sd of ? innovativeness')
hist(innov_mode,50);
ylabel('num of dyads')
xlabel('mode of ? innovativeness')
hist(explor_mean,50);
ylabel('num of dyads')
xlabel('mean of ? exploration')
hist(explor_std,50);
ylabel('num of dyads')
xlabel('sd of ? exploration')
hist(explor_mode,50);
ylabel('num of dyads')
xlabel('mode of ? exploration')

dyads_scores = [mf_dyads m_innov f_innov m_explor f_explor];
csvwrite('scores_by_dyads.csv',dyads_scores);

scores_by_member = [m_list innov_mode innov_std innov_mean explor_mode explor_std explor_mean];
csvwrite('scores_by_member.csv',scores_by_member);

mean(scores_by_member)
