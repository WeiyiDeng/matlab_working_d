% rng(1)
for t = 1:5
    for d = 1:10
        rng(1)
        sth = rand(1,1);
        disp(sth)
    end
end

% clearvars -EXCEPT rng