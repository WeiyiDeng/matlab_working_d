% scratchFolder = tempdir;
% numRows = 1e8;
% numColumns = 32;
% 
% filename = ['mmf' int2str(numRows) 'x' int2str(numColumns) '.dat'];
% filename = fullfile(scratchFolder, filename);
% 
% f = fopen(filename, 'w');
% for colNum = 1:numColumns
%     column = (1:numRows)' + colNum*1000000;
%     fwrite(f,column,'double');
% end
% fclose(f);
% 
% % % Prevent a memory-busting matrix from being created.
% % if numRows*numColumns*8 > 1e9
% %     error('Size possibly too big; are you sure you want to do this?')
% % end
% % 
% % mm = memmapfile(filename, 'Format', {'double', [numRows numColumns], 'm'});
% % m = mm.Data.m;
% % 
% % clear('m');
% 
% mm = memmapfile(filename, 'Format', {'double', [numRows 1], 'mj'}, 'Repeat', numColumns);
%            
% if ~isequal(mm.Data(8).mj, (1:numRows)' + 8*1000000)
%     error('The data was not read back in correctly!');
% end
% 
% summean = 0;
% for i = 1:numColumns
%     summean = summean + mean(mm.Data(i).mj);
% end

%%
% h = fopen('memmapex.raw','wb');
% for ii = 1:500
% fwrite(h,randn(60*256*256,1),'double');
% end
% fclose(h);
% 
% idx = randperm(30000);
% idx = sort(idx(1:10000));
%  
% rmean = zeros(256,256);
% m = memmapfile('memmapex.raw','Format',{'double',[256,256],'im'});
%  
% for ii = 1:length(idx)
% rmean = rmean + m.Data(idx(ii)).im;
% if mod(ii,100) == 0
% ii
% end
% end
%  
% rmean = rmean/length(idx);

%%
load('explor_contin_std.mat')
fileID = fopen('explor_contin_stdat.dat','w');
fwrite(fileID, explor_contin,'double');
fclose(fileID);

% mm = memmapfile('explor_contin_stdat.dat', 'Format', {'double', [47085403 1], 'col'},'Repeat', 6);
% % m = mm.Data(1).col;          % first col
% % info = whos('mm');
% m = mm.Data;
% 
% dtmean = [];
% for i = 1:8
%     dtmean = [dtmean mean(mm.Data(i).col)];
% end
% clearvars m mm



clear all

load('matpstd.mat');
IVs = matp(:,[6 7]);
clearvars matp

tic
beta_0 = [-5.6271    0.8375   27.4694    0.3870    0.1020];
beta_0 = [beta_0         -0.0509    0.0281   -0.2486    0.2183    0.1682   -0.1598    -0.0583    0.0421];  
beta_0 = [beta_0          0.0102    0.0140    0.0414   -0.0358    0.0538   -0.0406    0.0038   -0.0118];
beta_0 = [beta_0         -0.0112   -0.0003   -0.0245    0.0172    0.0076    0.0081   -0.0026    0.0057]; 
b = beta_0
b_explor = b([10:13 24:25 18:21 28:29])';

mm = memmapfile('explor_contin_stdat.dat', 'Format', {'double', [47085403 6], 'col'},'Repeat', 1);
% m = mm.Data(1).col;          % first col
% info = whos('mm');
m = mm.Data;

week_IV = 100*gampdf(IVs(:,2),b(2),b(3));          
week_IV(IVs(:,2)<1)=0;

explor_WD_multip = zeros(size(m.col));
for j = 1:6;
    explor_WD_multip(:,j) = m.col(:,j).*week_IV;
end

FV_explor = m.col*b_explor(1:6)+explor_WD_multip*b_explor(7:end);

clearvars m mm explor_WD_multip

toc

%%
clear all

load('matpstd.mat');
IVs = matp(:,[6 7]);
clearvars matp

tic
beta_0 = [-5.6271    0.8375   27.4694    0.3870    0.1020];
beta_0 = [beta_0         -0.0509    0.0281   -0.2486    0.2183    0.1682   -0.1598    -0.0583    0.0421];  
beta_0 = [beta_0          0.0102    0.0140    0.0414   -0.0358    0.0538   -0.0406    0.0038   -0.0118];
beta_0 = [beta_0         -0.0112   -0.0003   -0.0245    0.0172    0.0076    0.0081   -0.0026    0.0057]; 
b = beta_0
b_explor = b([10:13 24:25 18:21 28:29])';

week_IV = 100*gampdf(IVs(:,2),b(2),b(3));          
week_IV(IVs(:,2)<1)=0;

load('explor_contin_std.mat')

explor_WD_multip = zeros(size(explor_contin));
for j = 1:size(explor_contin,2);
    explor_WD_multip(:,j) = explor_contin(:,j).*week_IV;
end

FV_explor = explor_contin*b_explor(1:6)+explor_WD_multip*b_explor(7:end);
% FV_explor = [explor_contin explor_WD_multip]*b_explor;

clearvars m mm explor_WD_multip

toc