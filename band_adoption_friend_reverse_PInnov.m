clc
clear all

% adoptions = csvread('bandadoptions3.csv');               % note that adoptions has not subtracted 104
adoptions = csvread('bandadoptions_lenient_adopt.csv',1,0);    
adoptions = adoptions(:,1:3);
adopt_ones = ones(size(adoptions, 1),1);

% total number of weeks: length(105:527) = 423
T = 423
I = 8320                                                 % 165 members
J = 6046
old_bandt = 104
band_adoption = cell(T,1);

for t = 1:T
    [r c v] = find(adoptions(:,3)==(t+old_bandt));    
    adopt_ind = adoptions(r,:);
    band_adoption{t} = sparse(adopt_ind(:,1),adopt_ind(:,2),adopt_ones(r),I,J);
end

band_adopt_mat = sparse(adoptions(:,1),adoptions(:,2),adoptions(:,3)-old_bandt,I,J);                % important fix !!

friendlist = csvread('new_friendlist_8088.csv',1,0);

unique_member = unique(friendlist(:,1));
PI_member = zeros(size(unique_member));
for i = 1:length(unique_member)
    k_frineds_indices = find(friendlist(:,1)==unique_member(i));
    k_frineds_numbers = friendlist(k_frineds_indices,2);
    i_K_matrix = band_adopt_mat([unique_member(i); k_frineds_indices],:);
    count_i_J = sum(i_K_matrix(1,:)>0);
    Innov_i_no_SI = sum(i_K_matrix(1,:)>0 & i_K_matrix(1,:) == max(i_K_matrix));
    PI_member(i) = 1-Innov_i_no_SI/count_i_J;
end

PI_member = [unique_member PI_member];

save('PI_member.mat','PI_member') ;

