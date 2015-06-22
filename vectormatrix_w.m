A = [1 2 3;4 5 6]
repmat(A,2,3)                    % repeat the matrix twice in rows and 3 times in columns

A(:)                             % matrix vectorization

c = [1:12] 
reshape(c,3,4)                   % vector transformed into matrix
d = reshape(c,4,3)
d'

