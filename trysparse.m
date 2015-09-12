% load('member_dummies.mat');

% A = sparse([4:7],ones(1,4).*2,3,10,5);
% b = rand(5,1);
A = sparse([4:7],ones(1,4).*2,3,100000000,500);
b = rand(500,1);

% direct multiplication
tic
D = A*b;
toc

% compare with
% loop over columns (matlab stores sparse matrix in columns)
C = zeros(size(A,1),1);
tic
A = A.';
for i = 1:size(A,2)
    [r c v] = find(A(:,i));
    C(i) = sum(v.*b(r));
end
toc