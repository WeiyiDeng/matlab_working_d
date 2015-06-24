clc
clear

type = 1                            % choose 1 to input data, 2 to simulate, 3 to debug

if type == 1
    data = csvread('mod_data.csv');
    Xlocation = 1:24                 % columns of IVs     cannot be 0s
    ylocation = 29:32
    Ilocation = 25:28
    % b0 = zeros(1,18)
    b0 = zeros(1,20)                 % (J+K)*2
elseif type == 2
    % if simulation
    data = zeros(1000,3,1);          % IT*J*Kb
    Xlocation = 50                   % Xlocation is T now
    ylocation = []
    Ilocation = []
    % b0 = zeros(1,2)                  % (J-1+K)*2
    b0 = zeros(1,8)                  % (J+K)*2
    % b0(1) = 1
else                                 % type = 3
    % debug
    data = zeros(50,2,1);          % IT*J*K
    Xlocation = 5                   % Xlocation is T now
    ylocation = type                % do NOT modify
    Ilocation = []
    % b0 = zeros(1,8)                 % (J-1+K)*2
    b0 = zeros(1,6)                  % (J+K)*2
end

[betas, se, real_bs] = run_mnl(data, Xlocation, ylocation, Ilocation, b0)

disp(real_bs) 
disp(betas)
disp(se)