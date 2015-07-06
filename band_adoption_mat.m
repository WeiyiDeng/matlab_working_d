clc
clear

adoptions = csvread('bandadoptions.csv');
adopt_ones = ones(size(adoptions, 1),1);

% total number of weeks: length(53:527) = 475
T = 475
I = 194
J = 43666
band_adoption = cell(T,1);

for t = 1:T
    [r c v] = find(adoptions(:,3)==(t+52));    
    adopt_ind = adoptions(r,:);
    band_adoption{t} = sparse(adopt_ind(:,1),adopt_ind(:,2),adopt_ones(r),I,J);
end

% test
sum_a = 0;
for t = 1:T
    sum_a = sum_a + sum(sum([band_adoption{t}]));
end
sum_a

% cwd = pwd;
% cd(tempdir);
% pack                % pack can only be used in the command line ??
% cd(cwd)

for t = 1:T
    diffusion_jt(t,:) = sum(band_adoption{t},1);
end

% [ri,cj,vs] = find(diffusion_jt);
% [sm,sn] = size(diffusion_jt);
% diffusion_jt = sparse(ri,cj,vs,sm,sn);

cumu_diffusion = cumsum(diffusion_jt);


