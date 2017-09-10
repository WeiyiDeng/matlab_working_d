function [Q, gr, H] = obj_gmm(b_gmm,X,y,W0,n,k)

g = 1/n.*X'*(y-X*b_gmm);
Q = g'*W0*g;

end