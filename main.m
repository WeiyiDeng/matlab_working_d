clc
clear

type = 3                            % choose 1 to input data, 2 to simulate, 3 to debug

if type == 1
    data = csvread('mod_data.csv');
    Xlocation = 1:24                 % columns of IVs     cannot be 0s
    ylocation = 29:32
    Ilocation = 25:28
    b0 = zeros(1,18)
elseif type == 2
    % if simulation
    data = zeros(10000,5,0);          % IT*J*K
    Xlocation = 50                   % Xlocation is T now
    ylocation = []
    Ilocation = []
    b0 = zeros(1,8)
else                                 % type = 3
    % debug
    data = zeros(50,3,2);          % IT*J*K
    Xlocation = 5                   % Xlocation is T now
    ylocation = type
    Ilocation = []
    b0 = zeros(1,8)
end

[betas, se] = run_mnl(data, Xlocation, ylocation, Ilocation, b0)

disp(betas)
disp(se)