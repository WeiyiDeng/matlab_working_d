% load fisheriris
% X = meas(51:end,:);
% y = strcmp('versicolor',species(51:end));
% ds = mat2dataset(X);
% ds(1:5,:)
% ds.y = y;
% ds(1:5,:)
% modelspec = 'y ~ X1*X2*X3';
% mdl = fitglm(ds,modelspec,'Distribution','binomial')
% modelspec = 'y ~ X1+X2+X3';
% mdl = fitglm(ds,modelspec,'Distribution','binomial')
% modelspec = 'y ~ X1+X2+X3+X4';
% mdl = fitglm(ds,modelspec,'Distribution','binomial')

clc
clear

load('matp_friend_reverse.mat');
predict_trend = csvread('predict_trend_log4061_lenient.csv',1,0);

index = find(ismember(matp(:,3),predict_trend(:,1)));

predict_trend_mat = sparse(predict_trend(:,1),predict_trend(:,2),predict_trend(:,3),...
    max(predict_trend(:,1)),max(predict_trend(:,2)));

matp = matp(index,:);

trend_hat = zeros(size(matp,1),1);
for r = 1:length(trend_hat)
    trend_hat(r) = predict_trend_mat(matp(r,3),matp(r,4));
end
% matp(:,6) = trend_hat;

band_topics = csvread('band_count_topics10.csv',1,0);
band_topics_mat = sparse(band_topics(:,1),1,band_topics(:,2),max(band_topics(:,1)),1);
topics_count = zeros(size(matp,1),1);
for r = 1:length(topics_count)
    topics_count(r) = band_topics_mat(matp(r,3));
end

band_birth = csvread('introdate3.csv',1,0);
band_birth_mat = sparse(band_birth(:,1),1,band_birth(:,2),max(band_birth(:,1)),1);
band_age = zeros(size(matp,1),1);
for r = 1:length(band_age)
    band_age(r) = matp(r,4)-band_birth_mat(matp(r,3));
end
sum(band_age<0)

hist(band_topics(:,2))
% hist(band_birth(:,2))

% matp = [member_p friend_p band_p timeobs DV prob_adopt_week new_week_diff A_week_ijt];

load('dummies40.mat');
% size(dummies40)
dummies40 = dummies40(index,:);
dummies_prep = dummies40.*repmat(matp(:,8),1,4);
clearvars dummies40
dummy_prep = dummies_prep(:,1) | dummies_prep(:,2);         % 1-20

X = matp(:,6);
% X = matp(:,[6 7]);
y = matp(:,5);

clearvars matp 

load('N_prep_for_matp_jobs_organize.mat');
load('combi_similarities.mat')

N_prep_for_matp_jobs_organize = N_prep_for_matp_jobs_organize(index,:);

[row col val] = find(N_prep_for_matp_jobs_organize==1);
store_r = zeros(length(row),1);
for i = 1:length(row)-1
    if row(i+1) == row(i)+2
        store_r(i) = 1;
    end
end
row0 = row(store_r==1)+1;
col0 = col(store_r==1);
val0 = zeros(size(col0));

row_store = [row; row0];
col_store = [col; col0];
val_store = [val; val0];

% val_pdf = normpdf(val_store);
% tic
% N_prep_for_matp_jobs_organize=sparse(row_store,col_store,val_pdf);
% toc

None0s_X_N = [row_store col_store val_store];
X_N = N_prep_for_matp_jobs_organize;
S = combi_similarities;

clearvars N_prep_for_matp_jobs_organize combi_similarities

ds = mat2dataset(X);
ds.trend_hat = trend_hat;
ds.week_IV = dummy_prep;
ds.band_age = band_age;
ds.topics_count = topics_count;
ds.y = y;
clearvars -except ds

ds(1:5,:)
% modelspec = 'y ~ X1*X2*X3';
% mdl = fitglm(ds,modelspec,'Distribution','binomial')
modelspec = 'y ~ X1+trend_hat+week_IV+band_age+topics_count+week_IV*band_age+week_IV*topics_count';
mdl = fitglm(ds,modelspec,'Distribution','binomial')

