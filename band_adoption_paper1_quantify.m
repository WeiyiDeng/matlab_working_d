mean(week_IV)
std(week_IV)

ind = find(week_IV);
sth = week_IV(ind);
max(sth)
min(sth)
mean(sth)
median(sth)
std(sth)

% test = 0:0.1:10*std(sth);
test = 0:0.1:30*std(week_IV);
store = zeros(size(test));

for i = 1:length(test)
    FV = [mean(trend_hat)  test(i)  test(i).^2  mean(IV_N_S)]*b_basic;
    exp_util = exp(-(const+FV));
    prob=1./(1+exp_util);
    store(i) = prob;
end

plot(test,store)

% FV = [mean(trend_hat)  mean(week_IV)+std(week_IV)  (mean(week_IV)+std(week_IV)).^2  mean(IV_N_S)]*b_basic;
FV = [mean(trend_hat)  mean(week_IV)  mean(week_IV).^2   mean(pop)+100*std(pop)   mean(week_IV).*(mean(pop)+100*std(pop))      mean(week_IV).^2.*(mean(pop)+100*std(pop))    mean(IV_N_S)]*b_basic;
exp_util = exp(-(const+FV));
prob=1./(1+exp_util)

%%
test = 0:0.001:1*std(pop);
store = zeros(size(test));
for i = 1:length(test)
    FV = [mean(trend_hat)  mean(week_IV)  0.001.^2   test(i)   0.001.*test(i)      0.001.^2.*test(i)    mean(IV_N_S)]*b_basic;
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
    FV = [mean(trend_hat)  mean(week_IV)  mean(week_IV).^2   test(i)   mean(week_IV).*test(i)      mean(week_IV).^2.*test(i)    mean(IV_N_S)]*b_basic;
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
    FV = [mean(trend_hat)  mean(week_IV)  (mean(week_IV)+1*std(week_IV)).^2   test(i)   (mean(week_IV)+1*std(week_IV)).*test(i)      (mean(week_IV)+1*std(week_IV)).^2.*test(i)    mean(IV_N_S)]*b_basic;
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
    FV = [mean(trend_hat)  test(i)  test(i).^2   0.1712+std(pop)   test(i).*(0.1712+std(pop))      test(i).^2.*(0.1712+std(pop))    mean(IV_N_S)]*b_basic;
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
    FV = [mean(trend_hat)  test(i)  test(i).^2   0.1712-std(pop)   test(i).*(0.1712-std(pop))      test(i).^2.*(0.1712-std(pop))    mean(IV_N_S)]*b_basic;
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
    FV = [mean(trend_hat)  test(i)  test(i).^2   0.1712   test(i).*0.1712      test(i).^2.*0.1712    mean(IV_N_S)]*b_basic;
    exp_util = exp(-(const+FV));
    prob=1./(1+exp_util);
    store(i) = prob;
end
plot(test,store)
legend({'high pop','low pop','median pop'},'Location','northwest')
hold off