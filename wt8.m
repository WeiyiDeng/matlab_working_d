clc
clear
% social influence + external good + social influnce only lasts for certain
% time period + share of listenings + different beta for each user

%% preferences
J=3;       %number of alternatives
%w = 8;
%slot = 420
%T=w*slot;     %number of choices
T=100;
N=T*J;    %number of brands * number of choices
I=4;       % two users
% b=[1 1 3 14.8];
G = 5;     % external good utility range
%P = slot      % time period of social influence
P=T

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
features=randn(T,J,I);                                          %TXJx2
b = exp(randn(4,I))
U_external = G.*rand(T,I);

%% friendship matrix
friRan = rand(I);
index_fri = find(friRan<0.5);
Bfmat = zeros(I);
Bfmat(index_fri)=1;
Cfmat = tril(Bfmat,-1);                             % lower triangular matrix with diagonal elements =0
friends = Cfmat+Cfmat';

%% actual simulation
kIndex = 1;
k2Index = 1;
for t=1:T
    for i=1:I
        utility=exp(b(1,i)+b(2,i)*features(t,:,i)+b(3,i)*total_share{i}+b(4,i)*total_friends_share{i});    %1XJ
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
            choice=choice==test;                                        % 1*J
            total_choices{i}=total_choices{i}+choice;                         %1xJ
            total_share{i}= total_choices{i}./sum(total_choices{i});
        else
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