clc
clear

mat_export = csvread('mat_export.csv',1,0);
% load('dummy_mat.mat');

% % w: take too long to finish, need to run in parallel
% row = zeros(size(mat_export,1),1);
% col = zeros(size(mat_export,1),1);
% val = zeros(size(mat_export,1),1);
% for i = 1:size(mat_export,1)
%     row(i) = find(dummy_mat(:,1) == mat_export(i,1) & dummy_mat(:,2)== mat_export(i,3)...
%         & dummy_mat(:,3) == mat_export(i,4));
%     col(i) = mat_export(i,2);
%     val(i) = mat_export(i,7);
% end

innov = csvread('EAi3_lenient.csv');

innov_m = zeros(size(mat_export,1),1);
innov_f = zeros(size(mat_export,1),1);
for i = 1:size(mat_export,1)
    innov_m(i) = innov(mat_export(i,1),2);
    innov_f(i) = innov(mat_export(i,2),2);
end

innov_diff = innov_m-innov_f;            % abs?

% mat_export_innov = [mat_export innov_diff];

mat_export_innov_mf = [mat_export innov_m innov_f];

% csvwrite('mat_export_innov.csv',mat_export_innov);
csvwrite('mat_export_innov_mf.csv',mat_export_innov_mf);

%% remove duplicate double quotes in emEditor
dummy_innov_mat = csvread('dummy_SI_innov_mat.csv');
temp = dummy_innov_mat(end,:);
dummy_innov_mat = [temp; dummy_innov_mat];
dummy_innov_mat(end,:) = [];


