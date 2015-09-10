% load('member_dummies.mat');

A = sparse([4:7],ones(1,4).*2,3,10,5);
b = rand(5,1);

tic
D = A*b;
toc

C = zeros(size(A,2),1);
tic
A = A.';
for i = 1:size(A,2)
    [r c v] = find(A(:,i));
    C(i) = sum(v.*b(r));
end
toc