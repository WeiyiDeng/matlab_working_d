clc
clear

%% use of cell                       w: notice the use of { } !
A = cell(2,1)
A{1} = rand(5)
A{2} = rand(4)
A
% A{2} = {2 3 4}

%% preferences
J=3;       %number of alternatives
T=100;     %number of choices
N=T*J;    %number of brands * number of choices
I=2;       % two users
b=[1 1 0.1 .01];

%% variable initialization
% rng(0);
total_choices_of_friends = cell(I,1);
%total_choices{i} = cell(I,1);
total_choices = cell(I,1);
for i=1:I
    total_choices{i} = zeros(1,J);
    % total_choices_of_friends = zeros(1,J);
end
% features=randn(T,J,1);                                          %TXJx2
features=randn(T,J,I);                                          %TXJx2

%% friendship matrix
% friends=?

%% actual simulation
for t=1:T
    for i=1:I
        % utilty=exp(b(1)+b(2)*features(t,:,1)+b(3)*total_choices{i}+b(4)*total_choices_of_friends{1+I-i});   %1XJ
        utilty=exp(b(1)+b(2)*features(t,:,i)+b(3)*total_choices{i});
        utilty_other=sum(utilty,2);                                 %1x1        utility of all alternatives
        utilty_other=repmat(utilty_other,1,J);                      %1xJ
        prob=utilty./utilty_other;                                  %1xJ
        prob=cumsum(prob')';
        draw_for_choice=rand(1,1);
        draw_for_choice=repmat(draw_for_choice,1,J);                %1xJ
        choice=prob<draw_for_choice;
        choice=sum(choice,2)+1;                                     % 1*1
        choice=repmat(choice,1,J);                                  
        test=1:J;
        choice=choice==test;                                        % 1*J
        total_choices{i}=total_choices{i}+choice;                         %1xJ
        
       %% update total choices of friends
       % total_choices_of_friends{i}=?
       
    end
end

total_choices{1}
total_choices{2}