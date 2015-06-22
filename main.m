clc
clear

% data = csvread('mod_data.csv');
% Xlocation = 1:24                 % columns of IVs     cannot be 0s
% ylocation = 29:32
% Ilocation = 25:28
% b0 = zeros(1,18)

% if simulation
data = zeros(1000,5,2);          % IT*J*K
Xlocation = 10                   % Xlocation is T now
ylocation = []
Ilocation = []
b0 = zeros(1,12)

[betas, se] = run_mnl(data, Xlocation, ylocation, Ilocation, b0)

disp(betas)
disp(se)