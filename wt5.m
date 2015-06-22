clc
clear
% social influence + external good

%% preferences
J=3;       %number of alternatives
T=100;     %number of choices
N=T*J;    %number of brands * number of choices
I=4;       % two users
b=[1 1 0.1 .1];
G = 10     % external good utility range

%% variable initialization
rng(0);
total_choices_of_friends = cell(I,1);
total_choices = cell(I,1);
for i=1:I
    total_choices{i} = zeros(1,J);
    total_choices_of_friends{i} = zeros(1,J);
end
features=randn(T,J,I);                                          %TXJx2

U_external = G.*rand(T,I)

%% friendship matrix
friRan = rand(I)
index_fri = find(friRan<0.2)
Bfmat = zeros(I)
Bfmat(index_fri)=1
Cfmat = tril(Bfmat,-1)                             % lower triangular matrix with diagonal elements =0
friends = Cfmat+Cfmat'

%% actual simulation
for t=1:T
    for i=1:I
        utility=exp(b(1)+b(2)*features(t,:,i)+b(3)*total_choices{i}+b(4)*total_choices_of_friends{i});    %1XJ
        utility = [utility, exp(U_external(t,i))]
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
        else
        end
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

       % friends influence up to t-1
        choicemat = cell2mat(total_choices)
        for i = 1:I
            friendindex = repmat(friends(:,i),1,J)
            choicefriend = choicemat.*friendindex
            total_choices_of_friends{i} = sum(choicefriend,1)
        end
end

total_choices