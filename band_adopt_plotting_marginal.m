%% marginal plots !!
% WD set to 2
clear all

load('matpstd2.mat');
load('innov_contin_std2.mat');
load('explor_contin_std2.mat');

% sum(innov_contin(:,1)^2~=innov_contin(:,2))           % testing

mean_base_prob = mean(matp(:,6));                       % sd = 1
mean_wd = mean(matp(:,7))                              
std_wd = std(matp(:,7));
mean_innov_m = mean(innov_contin(:,1));                 % sd = 1
mean_innov_f = mean(innov_contin(:,3));                 % sd = 1
mean_explor_m = mean(explor_contin(:,1));               % sd = 1
mean_explor_f = mean(explor_contin(:,3));               % sd = 1

b = [-5.0951    0.19983    2.0459    0.3222    0.0627];
b = [b     -0.0102    0.0667   -0.1275    0.0555   -0.0005   -0.0652    -0.0225    0.0052];  
b = [b      0.0181   -0.0185    0.0131   -0.0125   -0.0088    0.0092   -0.0145    0.0017];
b = [b     -0.0160   -0.0046    0.0027    0.0005    0.0117    0.0008    0.0052   -0.0014]

const = b(1);

b_basic = b(4:5)';
b_innov = b([6:9 22:23 14:17 26:27])';
b_explor = b([10:13 24:25 18:21 28:29])';

% clear all
% load('innov_contin2.mat');
% load('explor_contin2.mat');

mean(innov_contin)
mean(explor_contin)
quantile(innov_contin(:,1),[0 0.05 0.1 0.25 0.3 0.5 0.6 0.75 0.9 0.95 1])
quantile(innov_contin(:,3),[0 0.05 0.1 0.25 0.3 0.5 0.6 0.75 0.9 0.95 1])
quantile(explor_contin(:,1),[0 0.05 0.1 0.25 0.3 0.5 0.6 0.75 0.9 0.95 1])
quantile(explor_contin(:,3),[0 0.05 0.1 0.25 0.3 0.5 0.6 0.75 0.9 0.95 1])

histogram(innov_contin(:,1),50)
histogram(innov_contin(:,3),50)
histogram(explor_contin(:,1),50)
histogram(explor_contin(:,3),50)
text(-1,120000,'14%')
text(2,120000,'86%')
line([1.9 1.9],[0 16*10^5], 'Color', 'r','Linestyle',':')
line([-0.84 -0.84],[0 16*10^5], 'Color', 'r','Linestyle','--')

%
q_eval = 0:0.01:1;
innov_m_q = quantile(innov_contin(:,1),q_eval);
innov_f_q = quantile(innov_contin(:,3),q_eval);
explor_m_q = quantile(explor_contin(:,1),q_eval);
explor_f_q = quantile(explor_contin(:,3),q_eval);

% WD_plug = mean_wd;
WD_plug = 2;
base_prob_plug = mean_base_prob      %+1;
innov_m_plug = mean_innov_m;
innov_f_plug = mean_innov_f;
explor_m_plug = mean_explor_m;
explor_f_plug = mean_explor_f;

prob = zeros(length(q_eval),1);
for k = 1:length(q_eval)
%     innov_m_plug = innov_m_q(k);
%     innov_f_plug = innov_f_q(k);                                             %<- to change
%     explor_m_plug = explor_m_q(k);
    explor_f_plug = explor_f_q(k);
    
    innov_plug = [innov_m_plug innov_m_plug^2 innov_f_plug innov_f_plug^2 innov_m_plug*innov_f_plug (innov_m_plug*innov_f_plug)^2];
    explor_plug = [explor_m_plug explor_m_plug^2 explor_f_plug explor_f_plug^2 explor_m_plug*explor_f_plug (explor_m_plug*explor_f_plug)^2];
    
    week_IV = 100*gampdf(WD_plug,exp(b(2)),exp(b(3)));
    week_IV(WD_plug<1)=0;
    
    innov_WD_multip = zeros(size(innov_plug));
    for i = 1:size(innov_plug,2);
        innov_WD_multip(:,i) = innov_plug(:,i).*week_IV;
    end
    
    explor_WD_multip = zeros(size(explor_plug));
    for j = 1:size(explor_plug,2);
        explor_WD_multip(:,j) = explor_plug(:,j).*week_IV;
    end
    
    FV = [base_prob_plug week_IV]*b_basic + innov_plug*b_innov(1:6)+innov_WD_multip*b_innov(7:end) + explor_plug*b_explor(1:6)+explor_WD_multip*b_explor(7:end);
    
    % % another way of calculating marginal effect of baseline probability
    % exp_util = exp(const+FV);
    % prob = exp_util/(1+exp_util)
    %
    % marginal_effect_based_prob = exp_util/(1+exp_util)^2*b(4)
    
    exp_util = exp(-(const+FV));         % this is now the utility of the external good
    prob(k)=1./(1+exp_util);
    % this is still the probability of choosing the product
    % pmat = [prob 1-prob];
    % pmat = pmat.*choice_dv;
    % [r c p] = find(pmat);                                             % I*1
    % ll_ur = sum(log(p))                 % ll of unrestricted model (the one that was estimated)
end


ax1 = axes;                               % generate two axes at same position
ax2 = axes('Position', get(ax1, 'Position'),'Color','none');

set(ax2,'YAxisLocation','right')          % move second axis to the right, remove x-ticks and labels
set(ax2,'XTick',[])
% set(ax2,'Xticklabel',[])

% histogram(ax1,innov_contin(:,1),30)
% histogram(innov_contin(:,1),50)
% histogram(ax1,innov_contin(:,3),30)                                          %<- to change
histogram(ax1,explor_contin(:,3),30)
% histogram(ax1,explor_contin(:,3),30)

hold on
% plot(innov_m_q,prob)

% plot(ax2,innov_f_q,prob,'r','linewidth',0.001)                               %<- to change
plot(ax2,explor_f_q,prob,'-r.','linewidth',0.001)
% plot(ax2,innov_m_q,prob,'r.')
hold off

x_point = quantile(explor_contin(:,3),[0.14 0.86]);                           %<- to change
% x_point = x_point +[0.06, -0.06];
% line([x_point(1) x_point(1)],[0 16*10^5], 'Color', 'r','Linestyle',':')
% line([x_point(2) x_point(2)],[0 16*10^5], 'Color', 'r','Linestyle','--')
line([x_point(1) x_point(1)],[0 0.014], 'Color', 'r','Linestyle','--')
line([x_point(2) x_point(2)],[0 0.014], 'Color', 'r','Linestyle','--')
text(x_point(1)+0.1,0.001,'14%')
text(x_point(2)+0.2,0.001,'86%')

ylabel(ax1,'hist of obs')
ylabel(ax2,'likelihood to adopt')
xlabel(ax1,'friend explor score')                                             %<- to change