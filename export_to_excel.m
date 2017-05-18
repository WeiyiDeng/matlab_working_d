%% testing
date_=clock;

LastName = {'Smith';'Johnson';'Williams';'Jones';'Brown'};                       % cell object here!
Age = [38;43;38;40;49];
Height = [71;69;64;67;64];
Weight = [176;163;131;133;119];
BloodPressure = [124 93; 109 77; 125 83; 117 75; 122 80];

T = table(Age,Height,Weight,BloodPressure,...
    'RowNames',LastName)

load('b.mat');
load('standard_error.mat');
load('t_stat.mat');
load('exit_flag.mat');

betaCoefficient = b';
StandardError = standard_error;
tStatistics = t_stat';
% Placeholder = {'1';'2';'3';'4';'5';'6';'7';'8';'9';'10';'11';'12';'13';...         % cell object here !
%     '14';'15';'16';'17';'18';'19';'20';'21';'22';'23';'24';'25';'26';'27';...
%     '28';'29'};
Placeholder = {'constant';'gamma_k';'gamma_?';'BP';'SI';'Innov_m';'Innov_m2';'Innov_f';'Innov_f2';...
    'Explor_m';'Explor_m2';'Explor_f';'Explor_f2';'Innov_m?Innov_f';'(Innov_m?Innov_f)2';...         % cell object here !
    'Explor_m?Explor_f';'(Explor_m?Explor_f)2';'Innov_m?SI';'Innov_m2?SI';'Innov_f?SI';'Innov_f2?SI';...
    'Explor_m?SI';'Explor_m2?SI';'Explor_f?SI';'Explor_f2?SI';'Innov_m?Innov_f?SI';'(Innov_m?Innov_f)2?SI';...
    'Explor_m?Explor_f?SI';'(Explor_m?Explor_f)2?SI'};

myTable = table(betaCoefficient, StandardError, tStatistics, 'RowNames', Placeholder)

mystr = strcat(num2str(date_(1)),num2str(date_(2)),num2str(date_(3)),num2str(date_(4)),num2str(date_(5)));
disp(mystr)

b_name = strcat('b_',mystr,'.mat');
SE_name = strcat('standard_error_',mystr,'.mat');
t_name = strcat('t_stat_',mystr,'.mat');
flag_name = strcat('exit_flag_',mystr,'.mat');

save(b_name,'b');
save(SE_name,'standard_error');
save(t_name,'t_stat');
save(flag_name,'exit_flag');

writetable(myTable,'myTable.xls','Sheet',1,'Range','B1','WriteRowNames',true)

% T = table(['M';'F';'M'],[45 45;41 32;40 34],...
%     {'NY';'CA';'MA'},[true;false;false])
% 
% writetable(T,'myData.xls','Sheet',2,'Range','B2:F6')

sth = Placeholder([1:5 29])                          % index some variable names from cell array

%% lenient sample export to xls
% stage 4
load('Results/mat_files/b_201733133')
load('Results/mat_files/t_stat_201733133')

betaCoefficient = b';
tStatistics= t_stat';

var_names = {'constant';'gamma_k';'gamma_*';'trend';'SI';'Innov_m';'Innov_m2';'Innov_f';'Innov_f2';...
    'Explor_m';'Explor_m2';'Explor_f';'Explor_f2';'Innov_m*Innov_f';'(Innov_m*Innov_f)2';...         % cell object here !
    'Explor_m*Explor_f';'(Explor_m*Explor_f)2';'Innov_m*SI';'Innov_m2*SI';'Innov_f*SI';'Innov_f2*SI';...
    'Explor_m*SI';'Explor_m2*SI';'Explor_f*SI';'Explor_f2*SI';'Innov_m*Innov_f*SI';'(Innov_m*Innov_f)2*SI';...
    'Explor_m*Explor_f*SI';'(Explor_m*Explor_f)2*SI';'new_baseline_prob'};

resultTable = table(betaCoefficient, tStatistics, 'RowNames', var_names)

writetable(resultTable,'resultTable_lenient.xls','Sheet',1,'Range','B1','WriteRowNames',true)

% stage 1
load('Results/mat_files/b_20173302355')
load('Results/mat_files/t_stat_20173302355')

betaCoefficient = b';
tStatistics= t_stat';

% var_names = {'constant';'gamma_k';'gamma_*';'trend';'SI';'Innov_m';'Innov_m2';'Innov_f';'Innov_f2';...
%     'Explor_m';'Explor_m2';'Explor_f';'Explor_f2';'Innov_m*Innov_f';'(Innov_m*Innov_f)2';...         % cell object here !
%     'Explor_m*Explor_f';'(Explor_m*Explor_f)2';'Innov_m*SI';'Innov_m2*SI';'Innov_f*SI';'Innov_f2*SI';...
%     'Explor_m*SI';'Explor_m2*SI';'Explor_f*SI';'Explor_f2*SI';'Innov_m*Innov_f*SI';'(Innov_m*Innov_f)2*SI';...
%     'Explor_m*Explor_f*SI';'(Explor_m*Explor_f)2*SI';'new_baseline_prob'};

resultTable = table(betaCoefficient, tStatistics)

writetable(resultTable,'resultTable_lenient.xls','Sheet',1,'Range','E1','WriteRowNames',false)

% stage 2
load('Results/mat_files/b_201733107')
load('Results/mat_files/t_stat_201733107')

betaCoefficient = b';
tStatistics= t_stat';

% var_names = {'constant';'gamma_k';'gamma_*';'trend';'SI';'Innov_m';'Innov_m2';'Innov_f';'Innov_f2';...
%     'Explor_m';'Explor_m2';'Explor_f';'Explor_f2';'Innov_m*Innov_f';'(Innov_m*Innov_f)2';...         % cell object here !
%     'Explor_m*Explor_f';'(Explor_m*Explor_f)2';'Innov_m*SI';'Innov_m2*SI';'Innov_f*SI';'Innov_f2*SI';...
%     'Explor_m*SI';'Explor_m2*SI';'Explor_f*SI';'Explor_f2*SI';'Innov_m*Innov_f*SI';'(Innov_m*Innov_f)2*SI';...
%     'Explor_m*Explor_f*SI';'(Explor_m*Explor_f)2*SI';'new_baseline_prob'};

resultTable = table(betaCoefficient, tStatistics)

writetable(resultTable,'resultTable_lenient.xls','Sheet',1,'Range','G1','WriteRowNames',false)

% stage 3
load('Results/mat_files/b_2017331022')
load('Results/mat_files/t_stat_2017331022')

betaCoefficient = b';
tStatistics= t_stat';

% var_names = {'constant';'gamma_k';'gamma_*';'trend';'SI';'Innov_m';'Innov_m2';'Innov_f';'Innov_f2';...
%     'Explor_m';'Explor_m2';'Explor_f';'Explor_f2';'Innov_m*Innov_f';'(Innov_m*Innov_f)2';...         % cell object here !
%     'Explor_m*Explor_f';'(Explor_m*Explor_f)2';'Innov_m*SI';'Innov_m2*SI';'Innov_f*SI';'Innov_f2*SI';...
%     'Explor_m*SI';'Explor_m2*SI';'Explor_f*SI';'Explor_f2*SI';'Innov_m*Innov_f*SI';'(Innov_m*Innov_f)2*SI';...
%     'Explor_m*Explor_f*SI';'(Explor_m*Explor_f)2*SI';'new_baseline_prob'};

resultTable = table(betaCoefficient, tStatistics)

writetable(resultTable,'resultTable_lenient.xls','Sheet',1,'Range','I1','WriteRowNames',false)

%% original sample export to xls
% stage 4
load('Results/mat_files/b_2017331646')
load('Results/mat_files/t_stat_2017331646')

betaCoefficient = b';
tStatistics= t_stat';

var_names = {'constant';'gamma_k';'gamma_*';'trend';'SI';'Innov_m';'Innov_m2';'Innov_f';'Innov_f2';...
    'Explor_m';'Explor_m2';'Explor_f';'Explor_f2';'Innov_m*Innov_f';'(Innov_m*Innov_f)2';...         % cell object here !
    'Explor_m*Explor_f';'(Explor_m*Explor_f)2';'Innov_m*SI';'Innov_m2*SI';'Innov_f*SI';'Innov_f2*SI';...
    'Explor_m*SI';'Explor_m2*SI';'Explor_f*SI';'Explor_f2*SI';'Innov_m*Innov_f*SI';'(Innov_m*Innov_f)2*SI';...
    'Explor_m*Explor_f*SI';'(Explor_m*Explor_f)2*SI';'new_baseline_prob'};

resultTable = table(betaCoefficient, tStatistics, 'RowNames', var_names)

writetable(resultTable,'resultTable_original.xls','Sheet',1,'Range','B1','WriteRowNames',true)

% stage 1
load('Results/mat_files/b_2017330186')
load('Results/mat_files/t_stat_2017330186')

betaCoefficient = b';
tStatistics= t_stat';

% var_names = {'constant';'gamma_k';'gamma_*';'trend';'SI';'Innov_m';'Innov_m2';'Innov_f';'Innov_f2';...
%     'Explor_m';'Explor_m2';'Explor_f';'Explor_f2';'Innov_m*Innov_f';'(Innov_m*Innov_f)2';...         % cell object here !
%     'Explor_m*Explor_f';'(Explor_m*Explor_f)2';'Innov_m*SI';'Innov_m2*SI';'Innov_f*SI';'Innov_f2*SI';...
%     'Explor_m*SI';'Explor_m2*SI';'Explor_f*SI';'Explor_f2*SI';'Innov_m*Innov_f*SI';'(Innov_m*Innov_f)2*SI';...
%     'Explor_m*Explor_f*SI';'(Explor_m*Explor_f)2*SI';'new_baseline_prob'};

resultTable = table(betaCoefficient, tStatistics)

writetable(resultTable,'resultTable_original.xls','Sheet',1,'Range','E1','WriteRowNames',false)

% stage 2
load('Results/mat_files/b_20173301847')
load('Results/mat_files/t_stat_20173301847')

betaCoefficient = b';
tStatistics= t_stat';

% var_names = {'constant';'gamma_k';'gamma_*';'trend';'SI';'Innov_m';'Innov_m2';'Innov_f';'Innov_f2';...
%     'Explor_m';'Explor_m2';'Explor_f';'Explor_f2';'Innov_m*Innov_f';'(Innov_m*Innov_f)2';...         % cell object here !
%     'Explor_m*Explor_f';'(Explor_m*Explor_f)2';'Innov_m*SI';'Innov_m2*SI';'Innov_f*SI';'Innov_f2*SI';...
%     'Explor_m*SI';'Explor_m2*SI';'Explor_f*SI';'Explor_f2*SI';'Innov_m*Innov_f*SI';'(Innov_m*Innov_f)2*SI';...
%     'Explor_m*Explor_f*SI';'(Explor_m*Explor_f)2*SI';'new_baseline_prob'};

resultTable = table(betaCoefficient, tStatistics)

writetable(resultTable,'resultTable_original.xls','Sheet',1,'Range','G1','WriteRowNames',false)

% stage 3
load('Results/mat_files/b_2017331020')
load('Results/mat_files/t_stat_2017331020')

betaCoefficient = b';
tStatistics= t_stat';

% var_names = {'constant';'gamma_k';'gamma_*';'trend';'SI';'Innov_m';'Innov_m2';'Innov_f';'Innov_f2';...
%     'Explor_m';'Explor_m2';'Explor_f';'Explor_f2';'Innov_m*Innov_f';'(Innov_m*Innov_f)2';...         % cell object here !
%     'Explor_m*Explor_f';'(Explor_m*Explor_f)2';'Innov_m*SI';'Innov_m2*SI';'Innov_f*SI';'Innov_f2*SI';...
%     'Explor_m*SI';'Explor_m2*SI';'Explor_f*SI';'Explor_f2*SI';'Innov_m*Innov_f*SI';'(Innov_m*Innov_f)2*SI';...
%     'Explor_m*Explor_f*SI';'(Explor_m*Explor_f)2*SI';'new_baseline_prob'};

resultTable = table(betaCoefficient, tStatistics)

writetable(resultTable,'resultTable_original.xls','Sheet',1,'Range','I1','WriteRowNames',false)


%% strict sample export to xls
% stage 4
load('Results/mat_files/b_2017331542')
load('Results/mat_files/t_stat_2017331542')

betaCoefficient = b';
tStatistics= t_stat';

var_names = {'constant';'gamma_k';'gamma_*';'trend';'SI';'Innov_m';'Innov_m2';'Innov_f';'Innov_f2';...
    'Explor_m';'Explor_m2';'Explor_f';'Explor_f2';'Innov_m*Innov_f';'(Innov_m*Innov_f)2';...         % cell object here !
    'Explor_m*Explor_f';'(Explor_m*Explor_f)2';'Innov_m*SI';'Innov_m2*SI';'Innov_f*SI';'Innov_f2*SI';...
    'Explor_m*SI';'Explor_m2*SI';'Explor_f*SI';'Explor_f2*SI';'Innov_m*Innov_f*SI';'(Innov_m*Innov_f)2*SI';...
    'Explor_m*Explor_f*SI';'(Explor_m*Explor_f)2*SI';'new_baseline_prob'};

resultTable = table(betaCoefficient, tStatistics, 'RowNames', var_names)

writetable(resultTable,'resultTable_strict.xls','Sheet',1,'Range','B1','WriteRowNames',true)

% stage 1
load('Results/mat_files/b_20173302223')
load('Results/mat_files/t_stat_20173302223')

betaCoefficient = b';
tStatistics= t_stat';

% var_names = {'constant';'gamma_k';'gamma_*';'trend';'SI';'Innov_m';'Innov_m2';'Innov_f';'Innov_f2';...
%     'Explor_m';'Explor_m2';'Explor_f';'Explor_f2';'Innov_m*Innov_f';'(Innov_m*Innov_f)2';...         % cell object here !
%     'Explor_m*Explor_f';'(Explor_m*Explor_f)2';'Innov_m*SI';'Innov_m2*SI';'Innov_f*SI';'Innov_f2*SI';...
%     'Explor_m*SI';'Explor_m2*SI';'Explor_f*SI';'Explor_f2*SI';'Innov_m*Innov_f*SI';'(Innov_m*Innov_f)2*SI';...
%     'Explor_m*Explor_f*SI';'(Explor_m*Explor_f)2*SI';'new_baseline_prob'};

resultTable = table(betaCoefficient, tStatistics)

writetable(resultTable,'resultTable_strict.xls','Sheet',1,'Range','E1','WriteRowNames',false)

% stage 2
load('Results/mat_files/b_20173302234')
load('Results/mat_files/t_stat_20173302234')

betaCoefficient = b';
tStatistics= t_stat';

% var_names = {'constant';'gamma_k';'gamma_*';'trend';'SI';'Innov_m';'Innov_m2';'Innov_f';'Innov_f2';...
%     'Explor_m';'Explor_m2';'Explor_f';'Explor_f2';'Innov_m*Innov_f';'(Innov_m*Innov_f)2';...         % cell object here !
%     'Explor_m*Explor_f';'(Explor_m*Explor_f)2';'Innov_m*SI';'Innov_m2*SI';'Innov_f*SI';'Innov_f2*SI';...
%     'Explor_m*SI';'Explor_m2*SI';'Explor_f*SI';'Explor_f2*SI';'Innov_m*Innov_f*SI';'(Innov_m*Innov_f)2*SI';...
%     'Explor_m*Explor_f*SI';'(Explor_m*Explor_f)2*SI';'new_baseline_prob'};

resultTable = table(betaCoefficient, tStatistics)

writetable(resultTable,'resultTable_strict.xls','Sheet',1,'Range','G1','WriteRowNames',false)

% stage 3
load('Results/mat_files/b_2017331151')
load('Results/mat_files/t_stat_2017331151')

betaCoefficient = b';
tStatistics= t_stat';

% var_names = {'constant';'gamma_k';'gamma_*';'trend';'SI';'Innov_m';'Innov_m2';'Innov_f';'Innov_f2';...
%     'Explor_m';'Explor_m2';'Explor_f';'Explor_f2';'Innov_m*Innov_f';'(Innov_m*Innov_f)2';...         % cell object here !
%     'Explor_m*Explor_f';'(Explor_m*Explor_f)2';'Innov_m*SI';'Innov_m2*SI';'Innov_f*SI';'Innov_f2*SI';...
%     'Explor_m*SI';'Explor_m2*SI';'Explor_f*SI';'Explor_f2*SI';'Innov_m*Innov_f*SI';'(Innov_m*Innov_f)2*SI';...
%     'Explor_m*Explor_f*SI';'(Explor_m*Explor_f)2*SI';'new_baseline_prob'};

resultTable = table(betaCoefficient, tStatistics)

writetable(resultTable,'resultTable_strict.xls','Sheet',1,'Range','I1','WriteRowNames',false)