load('matp_friend_reverse.mat');
predict_trend = csvread('predict_trend_log4061_lenient.csv',1,0);

predict_trend = csvread('predict_trend_log4061_lenient.csv',1,0);

index = find(ismember(matp(:,3),predict_trend(:,1)));

matp = matp(index,:);

sum(matp(:,7)>0)
size(matp)
ww = matp(matp(:,7)>0,7);
hist(ww)
dum1 = matp(:,7)<=10 & matp(:,7)>0;
dum2 = matp(:,7)<=20 & matp(:,7)>10;
dum3 = matp(:,7)<=30 & matp(:,7)>20;
dum4 = matp(:,7)<=40 & matp(:,7)>30;
dummies40 = [dum1 dum2 dum3 dum4];
dum1 = matp(:,7)<5;
dum2 = matp(:,7)<=10 & matp(:,7)>5;
dum3 = matp(:,7)<=15 & matp(:,7)>10;
dum4 = matp(:,7)<=30 & matp(:,7)>15;
dummies30 = [dum1 dum2 dum3 dum4];
save('dummies30f.mat','dummies30','-v7.3');
save('dummies40f.mat','dummies40','-v7.3');

load('dummies40f.mat');
% size(dummies40)
% dummies40 = dummies40(index,:);
dummies_prep = dummies40.*repmat(matp(:,8),1,4);
clearvars dummies40
dummy_prep = dummies_prep(:,1) | dummies_prep(:,2);         % 1-20

% mat_export = [matp(:,1:6) dummy_prep];
% csvwrite('mat_export.csv',mat_export)

csvwrite('dummy_prep_fix.csv',dummy_prep)
mat_export_fix = [matp(:,1:6) dummy_prep];
csvwrite('mat_export_fix.csv',mat_export_fix)
