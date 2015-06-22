clc
clear
N=9;    %number of bands (alternatives) * number of choices to make           w: total number of observations
J=3;       %number of alternatives
T=N/J;      %number of choices
b=[1 1 1];
rng(0);
features=randn(N,3);                                        %Nx3              w: features are variables here
utilty=exp(b(1)+b(2)*features(:,1)+b(3)*features(:,3));     %Nx1
utilty = reshape(utilty,T,J);                               %TxJ
utilty_other=sum(utilty,2);                                 %Tx1
utilty_other=repmat(utilty_other,1,J);                      %TxJ
prob=utilty./utilty_other;                                  %TxJ
prob=cumsum(prob')';
draw_for_choice=rand(T,1);
draw_for_choice=repmat(draw_for_choice,1,J);                %TxJ
choice=prob<draw_for_choice;
choice=sum(choice,2)+1;
choice=repmat(choice,1,J);
test=repmat(1:J,T,1);
choice=choice==test;