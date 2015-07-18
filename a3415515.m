clc
clear

% profile on


adoptions = csvread('bandadoptions.csv');
adopt_ones = ones(size(adoptions, 1),1);

% total number of weeks: length(105:527) = 423
T = 423
I = 194
J = 7600
old_bandt = 104
band_adoption = cell(T,1);

% profile viewer
% p = profile('info');
% 
% p.FunctionTable.ExecutedLines

% see http://nl.mathworks.com/help/matlab/ref/profile.html

