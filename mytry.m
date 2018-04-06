% % rng(1)
% for t = 1:4
% %     seed1 = 1; 
%     seed2 = 1;
%     for d = 1:2
% %         rng(100000+seed1)
%         for c = 1:3
%             rng(1000+seed2)
%             sth = rand(1,1);
%             seed2 = seed2+1;
%             disp(sth)
%         end
%     end
% end
% 
% % clearvars -EXCEPT rng

NB = 10;
T  = 52;
NM = 10;
beta_C = 0.3;
se_C = 10;

store = [];
for w=1:100
    sth = beta_C+se_C*randn(NB,T,NM);
    temp = mean(mean(mean(sth)));
    disp(temp)
    store = [store temp];
end

hist(store)
disp(mean(store))