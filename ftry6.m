% w: significance of innov member and friend ?
%% coefficients
clear all

% load('matpstd2.mat');
% load('innov_contin_std2.mat');
% load('explor_contin_std2.mat');
% 
% % sum(innov_contin(:,1)^2~=innov_contin(:,2))           % testing
% 
mean_base_prob = -0.2466;                       % sd = 1
% mean_wd = mean(matp(:,7))                              
% std_wd = std(matp(:,7));
% mean_innov_m = mean(innov_contin(:,1));                 % sd = 1
% mean_innov_f = mean(innov_contin(:,3));                 % sd = 1
% mean_explor_m = mean(explor_contin(:,1));               % sd = 1
% mean_explor_f = mean(explor_contin(:,3));               % sd = 1

% b = [-5.0951    0.19983    2.0459    0.3222    0.0627];
% b = [b     -0.0102    0.0667   -0.1275    0.0555   -0.0005   -0.0652    -0.0225    0.0052];  
% b = [b      0.0181   -0.0185    0.0131   -0.0125   -0.0088    0.0092   -0.0145    0.0017];
% b = [b     -0.0160   -0.0046    0.0027    0.0005    0.0117    0.0008    0.0052   -0.0014]

LB = [-5.0951    0.19983    2.0459    0.3222    0.0627];
LB = [LB     -0.0304    0.051   -0.1397    0.05073   -0.0005   -0.0652    -0.0225    0.0052];  
LB = [LB      0.005   -0.029    0.0069   -0.017   -0.0088    0.0092   -0.0145    0.0017];
LB = [LB     -0.028   -0.0074    0.0027    0.0005    0.0045    -0.001    0.0052   -0.0014]

UB = [-5.0951    0.19983    2.0459    0.3222    0.0627];
UB = [UB     0.01    0.0824   -0.1152    0.0603   -0.0005   -0.0652    -0.0225    0.0052];  
UB = [UB      0.0312   -0.0085    0.01934   -0.0079   -0.0088    0.0092   -0.0145    0.0017];
UB = [UB     -0.0039   -0.0018    0.0027    0.0005    0    0    0.0052   -0.0014]



% clear all
% load('innov_contin2.mat');
% load('explor_contin2.mat');

% mean(innov_contin)
% mean(explor_contin)
% quantile(innov_contin(:,1),[0 0.05 0.1 0.25 0.3 0.5 0.6 0.75 0.9 0.95 1])
% quantile(innov_contin(:,3),[0 0.05 0.1 0.25 0.3 0.5 0.6 0.75 0.9 0.95 1])
% quantile(explor_contin(:,1),[0 0.05 0.1 0.25 0.3 0.5 0.6 0.75 0.9 0.95 1])
% quantile(explor_contin(:,3),[0 0.05 0.1 0.25 0.3 0.5 0.6 0.75 0.9 0.95 1])
% quantile(innov_contin(:,1),[0.16 0.5 0.84])
% quantile(innov_contin(:,3),[0.16 0.5 0.84])
% quantile(explor_contin(:,1),[0.16 0.5 0.84])
% quantile(explor_contin(:,3),[0.16 0.5 0.84])

% hist(innov_contin(:,1),30)
% hist(innov_contin(:,3),30)
% hist(explor_contin(:,1),30)
% hist(explor_contin(:,3),30)
% line([1.9 1.9],[0 16*10^5], 'Color', 'r','Linestyle',':')
% line([-0.84 -0.84],[0 16*10^5], 'Color', 'r','Linestyle','--')

% WD_plug = mean_wd;
% WD_plug = 0;
base_prob_plug = mean_base_prob      %+1;
% innov_m_plug = mean_innov_m;
% innov_f_plug = mean_innov_f;
% explor_m_plug = mean_explor_m;
% explor_f_plug = mean_explor_f;
% innov_intv_m = [-0.8901   -0.2244    1.7034];     % 16% 50% 84% quantile
% innov_intv_f = [-0.9770   -0.1805    1.0119];     % 16% 50% 84% quantile
% explor_intv_m = [-0.8495   -0.4930    1.9075];    % 16% 50% 84% quantile
% explor_intv_f = [-0.7377   -0.2972    0.6773];    % 16% 50% 84% quantile
% innov_intv_m = quantile(innov_contin(:,1),[0.16 0.5 0.84])     % 16% 50% 84% quantile
% innov_intv_f = quantile(innov_contin(:,3),[0.16 0.5 0.84])     % 16% 50% 84% quantile
% explor_intv_m = quantile(explor_contin(:,1),[0.16 0.5 0.84])    % 16% 50% 84% quantile
% explor_intv_f = quantile(explor_contin(:,3),[0.16 0.5 0.84])    % 16% 50% 84% quantile
% % innov_intv_m = innov_intv_f

% % Get row numbers for each individual
% uni_m = unique(matp(:,1));
% uni_m_row = zeros(size(uni_m));
% for i = 1:length(uni_m)
%     ind_temp = find(matp(:,1)==uni_m(i));
%     uni_m_row(i) = ind_temp(1);
% end
% 
% uni_f = unique(matp(:,2));
% uni_f_row = zeros(size(uni_f));
% for i = 1:length(uni_f)
%     ind_temp = find(matp(:,2)==uni_f(i));
%     uni_f_row(i) = ind_temp(1);
% end
% 
% % test
% A = matp(uni_f_row,2);
% B = matp(uni_m_row,1);
% sum(A == uni_f)
% sum(B == uni_m)
%
% save('uni_m_row.mat','uni_m_row');
% save('uni_f_row.mat','uni_f_row');
% save('uni_m.mat','uni_m');
% save('uni_f.mat','uni_f');

% load('uni_m_row.mat');
% load('uni_f_row.mat');
% uni_innov_m = innov_contin(uni_m_row,1);
% uni_explor_m = explor_contin(uni_m_row,1);
% uni_innov_f = innov_contin(uni_f_row,3);
% uni_explor_f = explor_contin(uni_f_row,3);
% 
% save('uni_innov_m.mat','uni_innov_m');
% save('uni_explor_m.mat','uni_explor_m');
% save('uni_innov_f.mat','uni_innov_f');
% save('uni_explor_f.mat','uni_explor_f');

load('uni_innov_m.mat');
load('uni_explor_m.mat');
load('uni_innov_f.mat');
load('uni_explor_f.mat');

% innov_intv_m = quantile(uni_innov_m,[0.16 0.5 0.84 0.975])
% innov_intv_f = quantile(uni_innov_f,[0.16 0.5 0.84 0.975])
% explor_intv_m = quantile(uni_explor_m,[0.16 0.5 0.84 0.975])
% explor_intv_f = quantile(uni_explor_f,[0.16 0.5 0.84 0.975])
% innov_intv_m = quantile([uni_innov_m' uni_innov_f'],[0.16 0.5 0.84])
% innov_intv_f = innov_intv_m
% explor_intv_m = quantile([uni_explor_m' uni_explor_f'],[0.16 0.5 0.84])
% explor_intv_f = explor_intv_m

mincut_innov = quantile([uni_innov_m' uni_innov_f'],0.01);
maxcut_innov = quantile([uni_innov_m' uni_innov_f'],0.99);
mincut_explor = quantile([uni_explor_m' uni_explor_f'],0.01);
maxcut_explor = quantile([uni_explor_m' uni_explor_f'],0.99);
innov_intv_m = linspace(mincut_innov, maxcut_innov,33);
innov_intv_f = innov_intv_m;
explor_intv_m = linspace(mincut_explor, maxcut_explor,33);
explor_intv_f = explor_intv_m;
% innov_intv_m = quantile([uni_innov_m' uni_innov_f'],0.01:0.02:0.99);
% innov_intv_f = innov_intv_m;
% explor_intv_m = quantile([uni_explor_m' uni_explor_f'],0.01:0.02:0.99);
% explor_intv_f = explor_intv_m;

% innov_intv_m = -1.7694:0.1:3.5801;
% innov_intv_f = -1.7267:0.1:10.0343;
% explor_intv_m = -0.9172:0.1:1.9075;
% explor_intv_f = -0.9380:0.1:7.5189;

% innov_intv_m = sort(uni_innov_m);
% innov_intv_f = sort(uni_innov_f);
% explor_intv_m = sort(uni_innov_m);
% explor_intv_f = sort(uni_innov_f);

nbs = [LB; UB];

prob = zeros(size(nbs,1),length(innov_intv_m),length(innov_intv_f));      % To change!!

% WD_plug = -10:1:40;
WD_plug = 2;
for r = 1:size(nbs,1)
    b = nbs(r,:);
    
    const = b(1);
    
    b_basic = b(4:5)';
    b_innov = b([6:9 22:23 14:17 26:27])';
    b_explor = b([10:13 24:25 18:21 28:29])';
    
    for q = 1:length(innov_intv_m)                                                 % To change!!
        for g = 1:length(innov_intv_f)                                             % To change!!
            %     innov_m_plug = -0.2244;
            %     innov_f_plug = -0.1805;
            innov_m_plug = innov_intv_m(q);
            innov_f_plug = innov_intv_f(g);
            %         innov_m_plug = median([uni_innov_m' uni_innov_f']);
            %         innov_f_plug = median([uni_innov_m' uni_innov_f']);
            %     innov_f_plug = 8;                          % upside down shape of prob ?
            %             explor_m_plug = explor_intv_m(q);
            %             explor_f_plug = explor_intv_f(g);
            explor_m_plug = median([uni_explor_m' uni_explor_f']);
            explor_f_plug = median([uni_explor_m' uni_explor_f']);
            %         explor_m_plug = -0.4930;            % choose median score for other variables
            %         explor_f_plug = -0.2972;
            
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
            prob(r,q,g)=1./(1+exp_util);
            % this is still the probability of choosing the product
            % pmat = [prob 1-prob];
            % pmat = pmat.*choice_dv;
            % [r c p] = find(pmat);                                             % I*1
            % ll_ur = sum(log(p))                 % ll of unrestricted model (the one that was estimated)
            
        end
    end
end
LB_prob = prob(1,:,:);
UB_prob = prob(2,:,:);
sig_prob = double(squeeze(LB_prob >UB_prob));
% surf(innov_intv_m, innov_intv_f, sig_prob)
% % sub_prob = prob(2,:,:)-prob(1,:,:);
% % % surf(squeeze(prob(2,:,:)))
% % surf(squeeze(sub_prob))
% xlabel('innov_m');
% ylabel('innov_f');
% zlabel('sig_prob');

% plot(WD_plug, prob(:,1),'r')
% hold on 
% plot(WD_plug, prob(:,2),'k')
% plot(WD_plug, prob(:,3),'b')
% % plot(WD_plug, prob(:,4),'g')
% % plot(WD_plug, prob(:,5),'y')
% % legend('sender innov = L','sender innov = M','sender innov = H', 'sender innov = SH')
% % legend('sender innov = SL', 'sender innov = L','sender innov = M','sender innov = H', 'sender innov = SH')
% legend('sender explor = L','sender explor = M','sender explor = H')
% title('receiver explor score = high')
% xlabel('week difference')
% ylabel('likelihood')
% set(gca, 'YLim',[0 0.014])                    % change y axis
% hold off

%% coefficients
clearvars -EXCEPT sig_prob LB_prob UB_prob innov_intv_m innov_intv_f explor_intv_m explor_intv_f uni_explor_m uni_explor_f uni_innov_m uni_innov_f

% load('matpstd2.mat');
% load('innov_contin_std2.mat');
% load('explor_contin_std2.mat');
% 
% % sum(innov_contin(:,1)^2~=innov_contin(:,2))           % testing
% 
mean_base_prob = -0.2466;                       % sd = 1
% mean_wd = mean(matp(:,7))                              
% std_wd = std(matp(:,7));
% mean_innov_m = mean(innov_contin(:,1));                 % sd = 1
% mean_innov_f = mean(innov_contin(:,3));                 % sd = 1
% mean_explor_m = mean(explor_contin(:,1));               % sd = 1
% mean_explor_f = mean(explor_contin(:,3));               % sd = 1

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

% mean(innov_contin)
% mean(explor_contin)
% quantile(innov_contin(:,1),[0 0.05 0.1 0.25 0.3 0.5 0.6 0.75 0.9 0.95 1])
% quantile(innov_contin(:,3),[0 0.05 0.1 0.25 0.3 0.5 0.6 0.75 0.9 0.95 1])
% quantile(explor_contin(:,1),[0 0.05 0.1 0.25 0.3 0.5 0.6 0.75 0.9 0.95 1])
% quantile(explor_contin(:,3),[0 0.05 0.1 0.25 0.3 0.5 0.6 0.75 0.9 0.95 1])
% quantile(innov_contin(:,1),[0.16 0.5 0.84])
% quantile(innov_contin(:,3),[0.16 0.5 0.84])
% quantile(explor_contin(:,1),[0.16 0.5 0.84])
% quantile(explor_contin(:,3),[0.16 0.5 0.84])

% hist(innov_contin(:,1),30)
% hist(innov_contin(:,3),30)
% hist(explor_contin(:,1),30)
% hist(explor_contin(:,3),30)
% line([1.9 1.9],[0 16*10^5], 'Color', 'r','Linestyle',':')
% line([-0.84 -0.84],[0 16*10^5], 'Color', 'r','Linestyle','--')

% WD_plug = mean_wd;
% WD_plug = 0;
base_prob_plug = mean_base_prob      %+1;
% innov_m_plug = mean_innov_m;
% innov_f_plug = mean_innov_f;
% explor_m_plug = mean_explor_m;
% explor_f_plug = mean_explor_f;
% innov_intv_m = [-0.8901   -0.2244    1.7034];     % 16% 50% 84% quantile
% innov_intv_f = [-0.9770   -0.1805    1.0119];     % 16% 50% 84% quantile
% explor_intv_m = [-0.8495   -0.4930    1.9075];    % 16% 50% 84% quantile
% explor_intv_f = [-0.7377   -0.2972    0.6773];    % 16% 50% 84% quantile
% innov_intv_m = quantile(innov_contin(:,1),[0.16 0.5 0.84])     % 16% 50% 84% quantile
% innov_intv_f = quantile(innov_contin(:,3),[0.16 0.5 0.84])     % 16% 50% 84% quantile
% explor_intv_m = quantile(explor_contin(:,1),[0.16 0.5 0.84])    % 16% 50% 84% quantile
% explor_intv_f = quantile(explor_contin(:,3),[0.16 0.5 0.84])    % 16% 50% 84% quantile
% % innov_intv_m = innov_intv_f

% % Get row numbers for each individual
% uni_m = unique(matp(:,1));
% uni_m_row = zeros(size(uni_m));
% for i = 1:length(uni_m)
%     ind_temp = find(matp(:,1)==uni_m(i));
%     uni_m_row(i) = ind_temp(1);
% end
% 
% uni_f = unique(matp(:,2));
% uni_f_row = zeros(size(uni_f));
% for i = 1:length(uni_f)
%     ind_temp = find(matp(:,2)==uni_f(i));
%     uni_f_row(i) = ind_temp(1);
% end
% 
% % test
% A = matp(uni_f_row,2);
% B = matp(uni_m_row,1);
% sum(A == uni_f)
% sum(B == uni_m)
%
% save('uni_m_row.mat','uni_m_row');
% save('uni_f_row.mat','uni_f_row');
% save('uni_m.mat','uni_m');
% save('uni_f.mat','uni_f');

% load('uni_m_row.mat');
% load('uni_f_row.mat');
% uni_innov_m = innov_contin(uni_m_row,1);
% uni_explor_m = explor_contin(uni_m_row,1);
% uni_innov_f = innov_contin(uni_f_row,3);
% uni_explor_f = explor_contin(uni_f_row,3);
% 
% save('uni_innov_m.mat','uni_innov_m');
% save('uni_explor_m.mat','uni_explor_m');
% save('uni_innov_f.mat','uni_innov_f');
% save('uni_explor_f.mat','uni_explor_f');

% load('uni_innov_m.mat');
% load('uni_explor_m.mat');
% load('uni_innov_f.mat');
% load('uni_explor_f.mat');

% innov_intv_m = quantile(uni_innov_m,[0.16 0.5 0.84 0.975])
% innov_intv_f = quantile(uni_innov_f,[0.16 0.5 0.84 0.975])
% explor_intv_m = quantile(uni_explor_m,[0.16 0.5 0.84 0.975])
% explor_intv_f = quantile(uni_explor_f,[0.16 0.5 0.84 0.975])
% innov_intv_m = quantile([uni_innov_m' uni_innov_f'],[0.16 0.5 0.84])
% innov_intv_f = innov_intv_m
% explor_intv_m = quantile([uni_explor_m' uni_explor_f'],[0.16 0.5 0.84])
% explor_intv_f = explor_intv_m

% mincut_innov = quantile([uni_innov_m' uni_innov_f'],0.01);
% maxcut_innov = quantile([uni_innov_m' uni_innov_f'],0.99);
% mincut_explor = quantile([uni_explor_m' uni_explor_f'],0.01);
% maxcut_explor = quantile([uni_explor_m' uni_explor_f'],0.99);
% innov_intv_m = linspace(mincut_innov, maxcut_innov,50);
% innov_intv_f = innov_intv_m;
% explor_intv_m = linspace(mincut_explor, maxcut_explor,50);
% explor_intv_f = explor_intv_m;
% innov_intv_m = quantile([uni_innov_m' uni_innov_f'],0.01:0.02:0.99)
% innov_intv_f = innov_intv_m
% explor_intv_m = quantile([uni_explor_m' uni_explor_f'],0.01:0.02:0.99)
% explor_intv_f = explor_intv_m

% innov_intv_m = -1.7694:0.1:3.5801;
% innov_intv_f = -1.7267:0.1:10.0343;
% explor_intv_m = -0.9172:0.1:1.9075;
% explor_intv_f = -0.9380:0.1:7.5189;

% innov_intv_m = sort(uni_innov_m);
% innov_intv_f = sort(uni_innov_f);
% explor_intv_m = sort(uni_innov_m);
% explor_intv_f = sort(uni_innov_f);

% WD_plug = -10:1:40;
WD_plug = [0 2];

prob = zeros(length(WD_plug),length(innov_intv_m),length(innov_intv_f));      % To change!!
for q = 1:length(innov_intv_m)                                                 % To change!!
    for g = 1:length(innov_intv_f)                                             % To change!!
        %     innov_m_plug = -0.2244; 
        %     innov_f_plug = -0.1805;
        innov_m_plug = innov_intv_m(q);
        innov_f_plug = innov_intv_f(g);
%         innov_m_plug = median([uni_innov_m' uni_innov_f']);
%         innov_f_plug = median([uni_innov_m' uni_innov_f']);
        %     innov_f_plug = 8;                          % upside down shape of prob ?
%             explor_m_plug = explor_intv_m(q);
%             explor_f_plug = explor_intv_f(g);
            explor_m_plug = median([uni_explor_m' uni_explor_f']);
            explor_f_plug = median([uni_explor_m' uni_explor_f']);
%         explor_m_plug = -0.4930;            % choose median score for other variables
%         explor_f_plug = -0.2972;
        
        innov_plug = [innov_m_plug innov_m_plug^2 innov_f_plug innov_f_plug^2 innov_m_plug*innov_f_plug (innov_m_plug*innov_f_plug)^2];
        explor_plug = [explor_m_plug explor_m_plug^2 explor_f_plug explor_f_plug^2 explor_m_plug*explor_f_plug (explor_m_plug*explor_f_plug)^2];
        
        for k = 1:length(WD_plug)
            week_IV = 100*gampdf(WD_plug(k),exp(b(2)),exp(b(3)));
            week_IV(WD_plug(k)<1)=0;
            
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
            prob(k,q,g)=1./(1+exp_util);
            % this is still the probability of choosing the product
            % pmat = [prob 1-prob];
            % pmat = pmat.*choice_dv;
            % [r c p] = find(pmat);                                             % I*1
            % ll_ur = sum(log(p))                 % ll of unrestricted model (the one that was estimated)
        end
        
    end
end
sub_prob = prob(2,:,:)-prob(1,:,:);
% surf(squeeze(prob(2,:,:)))
h = surf(innov_intv_m, innov_intv_f, squeeze(sub_prob));
% h = surf(innov_intv_m, innov_intv_f, squeeze(prob(2,:,:)));
% h = surf(squeeze(prob(2,:,:)));
% sth = get(h,'ZData');
% set(h,'ZData',sth-0.005) 
xlabel('innov m')
ylabel('innov f')
zlabel('sig prob');
hold on
contour(innov_intv_m, innov_intv_f, sig_prob.*0.01)
% contour(sig_prob.*0.01)
% surf(innov_intv_m, innov_intv_f, sig_prob.*0.012)
% plot3(innov_intv_m, innov_intv_f, sig_prob,'.r','MarkerSize',1)
hold off

% plot(WD_plug, prob(:,1),'r')
% hold on 
% plot(WD_plug, prob(:,2),'k')
% plot(WD_plug, prob(:,3),'b')
% % plot(WD_plug, prob(:,4),'g')
% % plot(WD_plug, prob(:,5),'y')
% % legend('sender innov = L','sender innov = M','sender innov = H', 'sender innov = SH')
% % legend('sender innov = SL', 'sender innov = L','sender innov = M','sender innov = H', 'sender innov = SH')
% legend('sender explor = L','sender explor = M','sender explor = H')
% title('receiver explor score = high')
% xlabel('week difference')
% ylabel('likelihood')
% set(gca, 'YLim',[0 0.014])                    % change y axis
% hold off