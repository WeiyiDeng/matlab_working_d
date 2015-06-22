J = 50
I = 100


covmprep = eye(J).*0.5;
covm = repmat(covmprep,I,I);
covm2 = covm + (eye(J*I).*0.5);
L = chol(covm2,'lower');
sdndraw = randn(J*I,1);
bmv = zeros(J*I,1) + L*sdndraw;
bmv_mat = (reshape(bmv,J,I))';
exp(bmv_mat)

d = eye(500);
S = sparse(d)
[i,j,s] = find(S);