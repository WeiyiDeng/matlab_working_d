% test

% retry_BP

% x = 0:0.2:15;
% y = chi2pdf(x,4);
% plot(x,y)
% hold on
% y = chi2pdf(x,3);
% plot(x,y)
% y = chi2pdf(x,2);
% plot(x,y)
% y = chi2pdf(x,1);
% plot(x,y)
% hold off
% 
% x = 0:0.2:15;
% y = gampdf(x,2,0);
% plot(x,y)
% hold off
% 
% x = 0:0.2:30;
% for k = 0:10
%     for theta = 1:5
%         y = gampdf(x,k,theta);
%         plot(x,y)
%         hold on
%     end
% end
% hold off

gamma_k = 0.1351
gamma_theta = 2.6231
x = 0:0.2:52;
y = gampdf(x,exp(gamma_k),exp(gamma_theta));
plot(x,y)
hold on
gamma_k = 0.2252
gamma_theta = 2.4367
x = 0:0.2:52;
y = gampdf(x,exp(gamma_k),exp(gamma_theta));
plot(x,y)
hold on
gamma_k = 0.1727
gamma_theta = 2.5329
x = 0:0.2:52;
y = gampdf(x,exp(gamma_k),exp(gamma_theta));
plot(x,y)
hold on

x = 0:0.2:52;
y = gampdf(x,1.2212,7.736118);
plot(x,y)
hold off

%% generate simu data
clc
clear

N=200000;    %number of bands (alternatives) * number of choices to make           w: total number of observations
J=2;       %number of alternatives
T=N/J;      %number of choices
b=[1 0.5 0.7];

rng(0);

features=randn(T,3)-5;                                        %Nx3              w: features are variables here
% rand_index = randsample(N,100);
utilty=exp(b(1)+b(2)*features(:,1)+b(3)*features(:,3));     %Nx1
% utilty = reshape(utilty,T,J);                               %TxJ
% utilty_all=sum(utilty,2);                                 %Tx1
% utilty_all=repmat(utilty_all,1,J);                      %TxJ
prob=utilty./(utilty+1);                                  %TxJ
prob = [prob 1-prob];
prob=cumsum(prob')';
draw_for_choice=rand(T,1);
draw_for_choice=repmat(draw_for_choice,1,J);                %TxJ
choice=prob<draw_for_choice;
choice=sum(choice,2)+1;
choice=repmat(choice,1,J);
test=repmat(1:J,T,1);
choice=choice==test;

sum(choice,1)
max(prob)
min(prob)

test_matrix = [features(:,1) features(:,3) choice(:,1)];
save('test_matp.mat','test_matrix') ;
