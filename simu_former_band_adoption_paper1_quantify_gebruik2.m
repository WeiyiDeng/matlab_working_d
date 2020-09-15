clc
clear

rng(181102)

load('matp_friend_reverse.mat');
predict_trend = csvread('predict_trend_log4061_lenient.csv',1,0);

index = find(ismember(matp(:,3),predict_trend(:,1)));

matp = matp(index,:);

% load('dummies40.mat');
% % size(dummies40)
% dummies40 = dummies40(index,:);
% dummies_prep = dummies40.*repmat(matp(:,8),1,4);
% clearvars dummies40
% dummy_prep = dummies_prep(:,1) | dummies_prep(:,2);         % 1-20
% 
% % mat_export = [matp(:,1:6) dummy_prep];
% % csvwrite('mat_export.csv',mat_export)
% 
% prep_matp = matp(:,[1 3 4]);
% 
% % remove duplicte commas in EMeditor
% dummy_mat = csvread('dummy_SI_mat.csv');
% dum_mat_prep = dummy_mat(:,1:3);
% 
% [~,indx]=ismember(dum_mat_prep,prep_matp,'rows');
% 
% save('indx.mat','indx') ;
% save('dummy_mat.mat','dummy_mat') ;
load('indx.mat');
load('dummy_mat_fix.mat');

matp = matp(indx,:);

temp = matp(17574610:end,:);
matp(17574610:end,:) = [];
matp = [temp; matp];

% X = dummy_mat(:, 5);
DV_paper1 = dummy_mat(:,4);
dummy_agg_SI = dummy_mat(:, 6);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
band_pop = csvread('POP_counts_year_final.csv',1,0);
band_pop(:,2) = band_pop(:,2)-104;

band_pop_adjust_1 = band_pop;
temp_ind = find(band_pop_adjust_1(:,2)==423);
band_pop_adjust_1(temp_ind,:) = [];
band_pop_adjust_1(:,2) = band_pop_adjust_1(:,2)+1;

band_pop_mat = sparse(band_pop(:,1),band_pop(:,2),band_pop(:,4),max(band_pop(:,1)),max(band_pop(:,2)));
pop = zeros(size(matp,1),1);
for r = 1:length(pop)
    pop(r) = band_pop_mat(matp(r,3),matp(r,4));
end

band_pop_mat_adjust_1 = sparse(band_pop_adjust_1(:,1),band_pop_adjust_1(:,2),band_pop_adjust_1(:,4),max(band_pop_adjust_1(:,1)),max(band_pop_adjust_1(:,2)));
pop_adjust_1 = zeros(size(matp,1),1);
for r = 1:length(pop_adjust_1)
    pop_adjust_1(r) = band_pop_mat_adjust_1(matp(r,3),matp(r,4));
end

pop_adjust_1 = pop_adjust_1./1000;

predict_trend_mat = sparse(predict_trend(:,1),predict_trend(:,2),predict_trend(:,3),...
    max(predict_trend(:,1)),max(predict_trend(:,2)));

trend_hat = zeros(size(matp,1),1);
for r = 1:length(trend_hat)
    trend_hat(r) = predict_trend_mat(matp(r,3),matp(r,4));
end
matp(:,6) = trend_hat;

corr(trend_hat,pop)

band_topics = csvread('band_count_topics15.csv',1,0);
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

band_count_tracks = csvread('BAND_count_TRACKS.csv',1,0);
band_tracks_mat = sparse(band_count_tracks(:,1),1,band_count_tracks(:,2),max(band_count_tracks(:,1)),1);
tracks_count = zeros(size(matp,1),1);
for r = 1:length(tracks_count)
    tracks_count(r) = matp(r,4)-band_tracks_mat(matp(r,3));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('N_prep_for_matp_jobs_organize.mat');
load('combi_similarities.mat')

N_prep_for_matp_jobs_organize = N_prep_for_matp_jobs_organize(index,:);
N_prep_for_matp_jobs_organize = N_prep_for_matp_jobs_organize(indx,:);

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

matp_mjt = matp(:,[1 3 4]);

clearvars N_prep_for_matp_jobs_organize combi_similarities matp

%%
load('b_agg_dummy_pop3_full_fix.mat')
% beta_0 = b
% b = [-8.0906    0.2663   -0.0070    0.8637    -0.1342    0.0048   0.0391    2.1956   0.7717]
b = [-7.8942    0.2609   -0.0061    0.9530    -0.0451    0.0012   0.0445    2.1977   0.752];
% clearvars b

b_basic = b(2:7)';
const = b(1);

week_IV = dummy_agg_SI;

val_pdf = 100*normpdf(None0s_X_N(:,3),0,b(8));
X_N = sparse(None0s_X_N(:,1),None0s_X_N(:,2),val_pdf,17617085,6222);
IV_N_S = X_N*S.^exp(b(9));

pop = pop./1000;

% FV = [IVs(:,1)  trend_hat  week_IV  band_age  topics_count...
%     band_age.*week_IV  topics_count.*week_IV   band_age.*topics_count...
%     band_age.*topics_count.*week_IV    IV_N_S]*b_basic;          % with both trend and BP as controls
FV = [week_IV  week_IV.^2   pop   week_IV.*pop      week_IV.^2.*pop     IV_N_S]*b_basic;
exp_util = exp(-(const+FV));         % this is now the utility of the external good
prob=1./(1+exp_util);    

FV_rm = [pop   IV_N_S]*b_basic([3 6]);
exp_util_rm = exp(-(const+FV_rm));         % this is now the utility of the external good
prob_rm=1./(1+exp_util_rm); 


% mean(week_IV)
% std(week_IV)
% 
% ind = find(week_IV);
% sth = week_IV(ind);
% max(sth)
% min(sth)
% mean(sth)
% median(sth)
% std(sth)
% 
% % test = 0:0.1:10*std(sth);
% test = 0:0.1:30*std(week_IV);
% store = zeros(size(test));
% 
% for i = 1:length(test)
%     FV = [mean(trend_hat)  test(i)  test(i).^2  mean(IV_N_S)]*b_basic;
%     exp_util = exp(-(const+FV));
%     prob=1./(1+exp_util);
%     store(i) = prob;
% end
% 
% plot(test,store)

load('b_agg_dummy_pop3_full.mat')
b = [-7.8942    0.2609   -0.0061    0.9530    -0.0451    0.0012   0.0445    2.1977   0.752];
const = b(1);
b_basic = b(2:7)';
% FV = [mean(trend_hat)  mean(week_IV)+std(week_IV)  (mean(week_IV)+std(week_IV)).^2  mean(IV_N_S)]*b_basic;
FV = [mean(week_IV)  mean(week_IV).^2   mean(pop)+100*std(pop)   mean(week_IV).*(mean(pop)+100*std(pop))      mean(week_IV).^2.*(mean(pop)+100*std(pop))    mean(IV_N_S)]*b_basic;
% FV = [mean(trend_hat)  mean(week_IV)  mean(week_IV).^2   mean(pop)   mean(week_IV).*mean(pop)      mean(week_IV).^2.*mean(pop)    mean(IV_N_S)]*b_basic;
% FV = [median(trend_hat)  median(week_IV)  median(week_IV).^2   median(pop)+1*std(pop)   median(week_IV).*(median(pop)+1*std(pop))      median(week_IV).^2.*(median(pop)+1*std(pop))    median(IV_N_S)]*b_basic;
% FV = [median(trend_hat)  median(week_IV)  median(week_IV).^2   median(pop)+100*std(pop)   median(week_IV).*(median(pop)+100*std(pop))      median(week_IV).^2.*(median(pop)+100*std(pop))    median(IV_N_S)]*b_basic;
% FV = [mean(trend_hat)  mean(week_IV)  mean(week_IV).^2   mean(pop)-100*std(pop)   0    0   mean(IV_N_S)]*b_basic;
% FV = [mean(trend_hat)  mean(week_IV)  mean(week_IV).^2   mean(pop)   0      0    mean(IV_N_S)]*b_basic;
exp_util = exp(-(const+FV));
prob=1./(1+exp_util)

%%
test = 0:0.001:1*std(pop);
store = zeros(size(test));
for i = 1:length(test)
    FV = [mean(week_IV)  0.001.^2   test(i)   0.001.*test(i)      0.001.^2.*test(i)    mean(IV_N_S)]*b_basic;
    exp_util = exp(-(const+FV));
    prob=1./(1+exp_util);
    store(i) = prob;
end
plot(test,store)
% legend('high pop','Location','northwest')
hold on
test = 0:0.001:1*std(pop);
store = zeros(size(test));
for i = 1:length(test)
    FV = [mean(week_IV)  mean(week_IV).^2   test(i)   mean(week_IV).*test(i)      mean(week_IV).^2.*test(i)    mean(IV_N_S)]*b_basic;
    exp_util = exp(-(const+FV));
    prob=1./(1+exp_util);
    store(i) = prob;
end
plot(test,store)
% legend('low pop','Location','northwest')
hold on
test = 0:0.001:1*std(pop);
store = zeros(size(test));
for i = 1:length(test)
    FV = [mean(week_IV)  (mean(week_IV)+1*std(week_IV)).^2   test(i)   (mean(week_IV)+1*std(week_IV)).*test(i)      (mean(week_IV)+1*std(week_IV)).^2.*test(i)    mean(IV_N_S)]*b_basic;
    exp_util = exp(-(const+FV));
    prob=1./(1+exp_util);
    store(i) = prob;
end
plot(test,store)
legend({'low SI','median SI','high SI'},'Location','northwest')
hold off


%%
test = 0:0.1:15*std(week_IV);
store = zeros(size(test));
for i = 1:length(test)
    FV = [test(i)  test(i).^2   0.1712+std(pop)   test(i).*(0.1712+std(pop))      test(i).^2.*(0.1712+std(pop))    mean(IV_N_S)]*b_basic;
    exp_util = exp(-(const+FV));
    prob=1./(1+exp_util);
    store(i) = prob;
end
plot(test,store)
% legend('high pop','Location','northwest')
hold on
test = 0:0.1:15*std(week_IV);
store = zeros(size(test));
for i = 1:length(test)
    FV = [test(i)  test(i).^2   0.1712-std(pop)   test(i).*(0.1712-std(pop))      test(i).^2.*(0.1712-std(pop))    mean(IV_N_S)]*b_basic;
    exp_util = exp(-(const+FV));
    prob=1./(1+exp_util);
    store(i) = prob;
end
plot(test,store)
% legend('low pop','Location','northwest')
hold on
test = 0:0.1:15*std(week_IV);
store = zeros(size(test));
for i = 1:length(test)
    FV = [test(i)  test(i).^2   0.1712   test(i).*0.1712      test(i).^2.*0.1712    mean(IV_N_S)]*b_basic;
    exp_util = exp(-(const+FV));
    prob=1./(1+exp_util);
    store(i) = prob;
end
plot(test,store)
legend({'high pop','low pop','median pop'},'Location','northwest')
hold off

%%
test = 0:0.1:2.5*std(week_IV);
store = zeros(size(test));
for i = 1:length(test)
    FV = [test(i)  test(i).^2   0.1712+2*std(pop)   test(i).*(0.1712+std(pop))      test(i).^2.*(0.1712+std(pop))    mean(IV_N_S)]*b_basic;
    exp_util = exp(-(const+FV));
    prob=1./(1+exp_util);
    store(i) = prob;
end
plot(test(1:21),store(1:21))
% legend('high pop','Location','northwest')
hold on
test = 0:0.1:2.5*std(week_IV);
store = zeros(size(test));
for i = 1:length(test)
    FV = [test(i)  test(i).^2   0.1712-std(pop)   test(i).*(0.1712-std(pop))      test(i).^2.*(0.1712-std(pop))    mean(IV_N_S)]*b_basic;
    exp_util = exp(-(const+FV));
    prob=1./(1+exp_util);
    store(i) = prob;
end
plot(test(1:21),store(1:21))
% legend('low pop','Location','northwest')
hold on
test = 0:0.1:2.5*std(week_IV);
store = zeros(size(test));
for i = 1:length(test)
    FV = [test(i)  test(i).^2   0.1712   test(i).*0.1712      test(i).^2.*0.1712    mean(IV_N_S)]*b_basic;
    exp_util = exp(-(const+FV));
    prob=1./(1+exp_util);
    store(i) = prob;
end
plot(test(1:21),store(1:21))
legend({'high pop','low pop','median pop'},'Location','northwest')
hold off

%%
test = 0:0.1:2.5*std(week_IV);
store = zeros(size(test));
for i = 1:length(test)
    FV = [test(i)  test(i).^2   0.1712+2*std(pop)   test(i).*(0.1712+std(pop))      test(i).^2.*(0.1712+std(pop))    mean(IV_N_S)]*b_basic;
    exp_util = exp(-(const+FV));
    prob=1./(1+exp_util);
    store(i) = prob;
end
plot(test(1:21),store(1:21),'--','color','k')
% legend('high pop','Location','northwest')
hold on
test = 0:0.1:2.5*std(week_IV);
store = zeros(size(test));
for i = 1:length(test)
    FV = [test(i)  test(i).^2   0.1712   test(i).*0.1712      test(i).^2.*0.1712    mean(IV_N_S)]*b_basic;
    exp_util = exp(-(const+FV));
    prob=1./(1+exp_util);
    store(i) = prob;
end
plot(test(1:21),store(1:21),'-','color','k')
% legend('low pop','Location','northwest')
hold on
test = 0:0.1:2.5*std(week_IV);
store = zeros(size(test));
for i = 1:length(test)
    FV = [test(i)  test(i).^2   0.1712-std(pop)   test(i).*(0.1712-std(pop))      test(i).^2.*(0.1712-std(pop))    mean(IV_N_S)]*b_basic;
    exp_util = exp(-(const+FV));
    prob=1./(1+exp_util);
    store(i) = prob;
end
plot(test(1:21),store(1:21),'-.','color','k')
legend({'high pop','median pop','low pop'},'Location','northwest')
hold off
xlabel('Social influence level (1 unit is 1 std of social influence)') 
ylabel('Hazard rate')

%%
% b = [-7.8942    0.2609   -0.0061    0.9530    -0.0451    0.0012   0.0445    2.1977   0.752];
b = [-0.78942    0.2609   -0.0061    0.9530    -0.0451    0.0012   0.0445    2.1977   0.752];
% this will show concave curve!
const = b(1);
b_basic = b(2:7)';

% mean_SI = mean(week_IV(week_IV>0));
% median_SI = median(week_IV(week_IV>0));
% std_SI = std(week_IV(week_IV>0));

mean_SI = mean(week_IV);
median_SI = median(week_IV);
std_SI = std(week_IV);

med_SI = mean_SI;
med_SI = median_SI;
high_SI = mean_SI+2*std_SI;
low_SI = mean_SI-2*std_SI;

% med_SI = mean(week_IV(week_IV>0))
% high_SI = med_SI+std(week_IV)
% low_SI = med_SI-std(week_IV)

% high_SI = quantile(week_IV(week_IV>0),0.9)
% med_SI = quantile(week_IV(week_IV>0),0.5);
% low_SI = 0;

% med_SI = mean(week_IV(week_IV>0));
% high_SI = mean(week_IV(week_IV>0))+2*std(week_IV(week_IV>0));
% low_SI = mean(week_IV(week_IV>0))-2*std(week_IV(week_IV>0));
% high_SI = mean(week_IV(week_IV>0))+std(week_IV(week_IV>0));
% low_SI = mean(week_IV(week_IV>0))-std(week_IV(week_IV>0));

med_pop = mean(pop)
% med_pop = mean(pop(pop>0))
high_pop = med_pop + 2*std(pop)
low_pop = med_pop - 2*std(pop)
% med_pop = median(pop);
% high_pop = mean(pop) + 2*std(pop)
% low_pop = mean(pop) - 2*std(pop)
% high_pop = mean(pop) + std(pop);
% low_pop = mean(pop) - std(pop);

% test = [low_SI med_SI high_SI];
test = [0 1 2];
pop_test = [low_pop med_pop high_pop];
store_u = zeros(length(test),length(pop_test));
store_p = zeros(length(test),length(pop_test));
for i = 1:length(test)
    for j = 1:length(pop_test)
        FV = [test(i)  test(i).^2   pop_test(j)   test(i).*pop_test(j)      test(i).^2.*pop_test(j)    mean(IV_N_S)]*b_basic;
        exp_util = exp(-(const+FV));
        prob=1./(1+exp_util);
        store_p(i,j) = prob;
        store_u(i,j) = [test(i)  test(i).^2   pop_test(j)   test(i).*pop_test(j)      test(i).^2.*pop_test(j)]*b_basic(1:5);
    end
end
store_p
store_u


%%
% med_SI = mean(week_IV(week_IV>0))
% high_SI = med_SI+std(week_IV)
% low_SI = med_SI-std(week_IV)

med_pop = mean(pop)
% med_pop = mean(pop(pop>0))
high_pop = med_pop + 2*std(pop)
low_pop = med_pop - 2*std(pop)

test = 0:0.5:8;
store = zeros(size(test));
store_u = zeros(size(test));
for i = 1:length(test)
    FV = [test(i)  test(i).^2   high_pop   test(i).*high_pop      test(i).^2.*high_pop    mean(IV_N_S)]*b_basic;
    exp_util = exp(-(const+FV));
    prob=1./(1+exp_util);
    store(i) = prob;
    store_u(i) = [test(i)  test(i).^2   high_pop   test(i).*high_pop      test(i).^2.*high_pop]*b_basic(1:5);
end
plot(test,store.*0.001,'--','color','k')
% plot(test,store_u,'--','color','k')
% legend('high pop','Location','northwest')
hold on
test = 0:0.5:8;
store = zeros(size(test));
store_u = zeros(size(test));
for i = 1:length(test)
    FV = [test(i)  test(i).^2   med_pop   test(i).*med_pop      test(i).^2.*med_pop    mean(IV_N_S)]*b_basic;
    exp_util = exp(-(const+FV));
    prob=1./(1+exp_util);
    store(i) = prob;
    store_u(i) = [test(i)  test(i).^2   med_pop   test(i).*med_pop      test(i).^2.*med_pop]*b_basic(1:5);
end
plot(test,store.*0.001,'-','color','k')
% plot(test,store_u,'-','color','k')
% legend('low pop','Location','northwest')
hold on
test = 0:0.5:8;
store = zeros(size(test));
store_u = zeros(size(test));
for i = 1:length(test)
    FV = [test(i)  test(i).^2   low_pop   test(i).*low_pop      test(i).^2.*low_pop    mean(IV_N_S)]*b_basic;
    exp_util = exp(-(const+FV));
    prob=1./(1+exp_util);
    store(i) = prob;
    store_u(i) = [test(i)  test(i).^2   low_pop   test(i).*low_pop      test(i).^2.*low_pop]*b_basic(1:5);
end
plot(test,store.*0.001,'-.','color','k')
% plot(test,store_u,'-.','color','k')
legend({'high pop','median pop','low pop'},'Location','northwest')
hold off
xlabel('Social influence (number of friends sampling a new band b in a week') 
ylabel('marginal contribution of social influence')
% xticks([0 10 20])
% xticklabels({'low','median','high'})

%%
med_pop = mean(pop)
% med_pop = median(pop)
% med_pop = mean(pop(pop>0))
test = 0:1:20;
store = zeros(size(test));
store_u = zeros(size(test));
for i = 1:length(test)
    FV = [test(i)  test(i).^2   med_pop   test(i).*med_pop      test(i).^2.*med_pop    mean(IV_N_S)]*b_basic;
    exp_util = exp(-(const+FV));
    prob=1./(1+exp_util);
    store(i) = prob;
    store_u(i) = [test(i)  test(i).^2   med_pop   test(i).*med_pop      test(i).^2.*med_pop]*b_basic(1:5);
end
figure
plot(test,store,'-','color','k')
% figure
% plot(test,store_u,'-','color','k')
xlabel('Social influence (number of friends sampling a new band b in a week)') 
ylabel('Hazard rate')

