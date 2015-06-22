clc
clear
% social influence + external good + social influnce only lasts for certain
% time period + squared choices + trend of whole market share + dying band
% + innovators who get tried fast 
% + different betas for each user and each band ?
% + way to draw from truncated normal +way to draw from multivariate normal
%% preferences
J=50;       %number of alternatives
%w = 8;
%slot = 420
%T=w*slot;     %number of choices
T=420;
N=T*J;    %number of brands * number of choices
I=10;       % two users
% b=[1 1 .1 .1 -0.001 1];
G = 5;     % external good utility range
%P = slot      % time period of social influence
P=420
Ps = 420        % time period to replace dying bands

%% variable initialization
rng(0);
total_choices_of_friends = cell(I,1);
total_choices = cell(I,1);
total_choices_archive = cell(I,1);
total_share = cell(I,1);
total_friends_share = cell(I,1);
for i=1:I
    total_choices{i} = zeros(1,J);
    total_choices_of_friends{i} = zeros(1,J);
    total_choices_archive{i} = zeros(P,J);
    total_share{i} = zeros(1,J);
    total_friends_share{i} = zeros(1,J);
end
ichoice = zeros(I,J);
trend = zeros(T,J);                                                             % out of memory why store T?
total_cshare_archive = zeros(Ps,J);

U_external = G.*rand(T,I);

features=randn(T,J,I);                                          %TXJx2              % out of memory T?

% b = exp(randn(4,I))                      % different betas for each user
% b = exp(rand(4,J))
h = [1 1 1 0.5 1 1]
beta1 = exp(randn(I,J)).*h(1);
beta2 = exp(randn(I,J)).*h(2);
beta3 = exp(randn(I,J)).*h(3);
beta4 = exp(randn(I,J)).*h(4);
beta6 = exp(randn(I,J)).*h(6);
% beta6 = exp(1+ 2.*randn(I,J))


%% draw beta5 from truncated normal 
F_a = normcdf(-3);
F_b = normcdf(3);
mu = rand(I,J);
mu_trun = (1-mu).*F_a + mu.*F_b;
btrun_norm = norminv(mu_trun);
beta5 = -exp(btrun_norm).*h(5);

% generate innovators with higher beta5
innovb = 10;
innovRan = rand(I,1);
index_innov = find(innovRan<0.1);
innovator = zeros(I,1);
innovator(index_innov) = 1;
beta5 = beta5.*repmat(1-innovator,1,J) + beta5.*repmat(innovator,1,J).*innovb;

%% trial with multivariate normal draws with covariance matrix
% covmprep = eye(J).*0.5;
% covm = repmat(covmprep,I,I);
% covm2 = covm + (eye(J*I).*0.5);
% L = chol(covm2,'lower');
% sdndraw = randn(J*I,1);
% bmv = zeros(J*I,1) + L*sdndraw;
% bmv_mat = (reshape(bmv,J,I))';
% exp(bmv_mat)

%% friendship matrix
friRan = rand(I);
index_fri = find(friRan<0.2);
Bfmat = zeros(I);
Bfmat(index_fri)=1;
Cfmat = tril(Bfmat,-1);                             % lower triangular matrix with diagonal elements =0
friends = Cfmat+Cfmat';

%% actual simulation
kIndex = 1;
k2Index = 1;
k3Index = 1;
old_mktshare = zeros(1,J);
for t=1:T
    for i=1:I
        % utility=exp(b(1,i)+b(2,i)*features(t,:,i)+b(3,i)*total_share{i}+b(4,i)*total_friends_share{i});    %1XJ
        % utility=exp(b(1,:)+b(2,:).*features(t,:,i)+b(3,:).*total_share{i}+b(4,:).*total_friends_share{i});    %1XJ
        % utility=exp(b(1)+b(2)*features(t,:,i)+b(3)*total_choices{i}+b(4)*total_choices_of_friends{i}+b(5)*(total_choices{i}).^2 +b(6)*trend(t,:));    %1XJ
        utility=exp(beta1(i,:)+beta2(i,:).*features(t,:,i)+beta3(i,:).*total_choices{i}+beta4(i,:).*total_choices_of_friends{i}+beta5(i,:).*(total_choices{i}).^2 +beta6(i,:).*trend(t,:));    %1XJ
        utility = [utility, exp(U_external(t,i))];
        utility_other=sum(utility,2);                                 %1x1        utility of all alternatives
        utility_other=repmat(utility_other,1,(J+1));                      %1xJ
        prob=utility./utility_other;                                  %1xJ
        prob=cumsum(prob')';
        draw_for_choice=rand(1,1);
        draw_for_choice=repmat(draw_for_choice,1,(J+1));                %1xJ
        choice=prob<draw_for_choice;
        choice=sum(choice,2)+1;                                     % 1*1
        if choice <= J
            choice=repmat(choice,1,J);
            test=1:J;
            ichoice(i,:)=choice==test;                                        % 1*J
            total_choices{i}=total_choices{i}+ichoice(i,:);                         %1xJ
            total_share{i}= total_choices{i}./sum(total_choices{i});
        else
            ichoice(i,:) = zeros(1,J);                              % beware change
        end
        % store histry of past choices of each time point
        if t <= P
            total_choices_archive{i}(t,:) = total_choices{i};
        else total_choices_archive{i}(t-kIndex*P,:) = total_choices{i};
        end
    end
    
    if t == (kIndex+1)*P
        kIndex = kIndex+1;
    else
    end
    
    %% fashion trend: whole market share of t minus share of t-1
    total_choicet = sum(ichoice,1);
    new_mktshare = total_choicet./sum(total_choicet);
    trend((t+1),:) = new_mktshare - old_mktshare;
    old_mktshare = new_mktshare;
    
    %% change bands in the choice set
    if t <= Ps
        total_cshare_archive(t,:) = total_choicet;
    else total_cshare_archive(t-k3Index*Ps,:) = total_choicet;
    end
    
    if t == (k3Index+1)*Ps
        k3Index = k3Index+1;
        share_in_ps = sum(total_cshare_archive,1)./sum(sum(total_cshare_archive));     % replace dying bands
        dyband = share_in_ps<0.01;
        choicematp2 = cell2mat(total_choices);
        dyprep = repmat(1-dyband,I,1);
        killband = choicematp2.*dyprep;
        total_choices = num2cell(killband,2);
        dyindex = find(dyband);
        for k= 1:length(dyindex)
            beta1(:,k) = exp(randn(I,1)).*h(1);
            beta2(:,k) = exp(randn(I,1)).*h(2);
            beta3(:,k) = exp(randn(I,1)).*h(3);
            beta4(:,k) = exp(randn(I,1)).*h(4);
            beta6(:,k) = exp(randn(I,1)).*h(6);
            mu = rand(I,1);
            mu_trun = (1-mu).*F_a + mu.*F_b;
            btrun_norm = norminv(mu_trun);
            beta5(:,k) = -exp(btrun_norm).*h(5);
            beta5(:,k) = beta5(:,k).*(1-innovator) + beta5(:,k).*innovator.*innovb;
        end
    else
    end
    
%     if t > Ps
%         share_in_ps = sum(total_cshare_archive,1)./sum(sum(total_cshare_archive))
%         dyband = share_in_ps<0.1
%         choicematp2 = cell2mat(total_choices)
%         dyprep = repmat(1-dyband,I,1)
%         killband = choicematp2.*dyprep
%         total_choices = num2cell(killband,2)
%     else
%     end
    
    %% update total choices of friends
    % simplified version
    %        choicemat = cell2mat(total_choices)
    %        friprep = repmat(friends(:,i),1,J)
    %        cfmultip = choicemat.*friprep
    %        total_choices_of_friends{i}= sum(cfmultip,1)
    
    % complex version (friends influence in real time)
    %        choicemat = cell2mat(total_choices)
    %        for i=1:I
    %            friendi = find(friends(:,i))
    %            if length(friendi) < 1,
    %            else
    %                for j = 1:length(friendi)
    %                    friprep = repmat(friends(:,friendi(j)),1,J)
    %                    cfmultip = choicemat.*friprep
    %                    total_choices_of_friends{friendi(j)}= sum(cfmultip,1)
    %                end
    %            end
    %        end
    
    %% friends influence up to t-1 (influence last for period P)
    
    % cellfun(@(x) x(10,1), total_choices_archive, 'UniformOutput', false)
    % total_choices_archive{:}(10,1) is not supported
    if t == ((k3Index-1)+1)*Ps && k3Index ~= 1
        for i = 1:I
            total_choices_archive{i}(:,dyindex) = zeros(P,length(dyindex));
        end
    end
         
    if t > P
        fcindex = t-k2Index*P+1;                    % former choice index
        former_choices = cellfun(@(x) x(fcindex,:), total_choices_archive, 'UniformOutput', false);
        former_choices = cell2mat(former_choices);
    else
        former_choices = zeros(I,J);
    end
    
    choicematp = cell2mat(total_choices);           % changed
    choicemat = choicematp-former_choices;           % calculate choices within past time period P
    
    for i = 1:I
        friendindex = repmat(friends(:,i),1,J);
        choicefriend = choicemat.*friendindex;
        total_choices_of_friends{i} = sum(choicefriend,1);
        total_friends_share{i} = total_choices_of_friends{i}./sum(total_choices_of_friends{i});
    end
    
    if t == (k2Index+1)*P-1
        k2Index = k2Index+1;
    else
    end
    
end

total_choices