% rng(1)
for t = 1:4
%     seed1 = 1; 
    seed2 = 1;
    for d = 1:2
%         rng(100000+seed1)
        for c = 1:3
            rng(1000+seed2)
            sth = rand(1,1);
            seed2 = seed2+1;
            disp(sth)
        end
    end
end

% clearvars -EXCEPT rng