clear all

load('matpstd2.mat');
load('uni_m_row.mat');
load('uni_f_row.mat');

A = sort(uni_m_row);
B = sort(uni_f_row);

overlap_bands = zeros(length(B),1);
for i = 1:length(B)-1
    tempvec = unique(matp(B(i):B(i+1)-1,3));
    overlap_bands(i) = length(tempvec);
end
overlap_bands(i+1) = length(unique(matp(B(i):end,3)));

mean(overlap_bands)