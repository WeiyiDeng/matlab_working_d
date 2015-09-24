
A = rand(100,1);
B = rand(100,1);
corr(A,B)
vecA = [];
vecB = [];
for i = 1:length(A)
    rtimes = round(rand(1)*1000000);
    vecA = [vecA; repmat(A(i),rtimes,1)];
    vecB = [vecB; repmat(B(i),rtimes,1)];
end
corr(vecA,vecB)

innov = innov(:,2);

list2 = []
for i = 2:size(innov_m,1)
    if innov_m(i,1)~=innov_m(i-1,1)
        list2 = [list2 i-1];
    else
    end
end
list2 = [list2 size(innov_m,1)];    

member_namelist = matp(list2,1);
small_in = innov(member_namelist);
small_exp = explor(member_namelist);
member_rowlist = list2(1);
corr(small_in,small_exp)
for i = 2:length(list2)
    member_rowlist = [member_rowlist (list2(i) - list2(i-1))];
end
try_in = []
try_exp = []
for i = 1:length(member_namelist)
    try_in = [try_in; repmat(small_in(i),member_rowlist(i),1)];
    try_exp = [try_exp; repmat(small_exp(i),member_rowlist(i),1)];
end
corr(try_in,try_exp)
save('member_rowlist.mat','member_rowlist')

load('member_rowlist.mat')
rng(1)
rv1 = rand(158,1);
rv2 = round(rand(158,1).*200);
corr(rv1,rv2)
try_rv1 = []
try_rv2 = []
for i = 1:length(member_namelist)
    try_rv1 = [try_rv1; repmat(rv1(i),member_rowlist(i),1)];
    try_rv2 = [try_rv2; repmat(rv2(i),member_rowlist(i),1)];
end
corr(try_rv1,try_rv2)

%%
% gamma function is differentiatable to x (may not wrt k and theta)
delta = eps
(gampdf(1,1+delta,2)-gampdf(1,1,2))/delta
delta = -delta
(gampdf(1,1+delta,2)-gampdf(1,1,2))/delta

delta = 10*eps
(gampdf(10,5+delta,6)-gampdf(10,5,6))/delta

delta = -10*eps
(gampdf(10,5+delta,6)-gampdf(10,5,6))/delta

((2+delta)^2-2^2)/delta

delta = 10*eps
((2+delta)^2-2^2)/delta

