function [Q, gr, H] = obj_gmm_IV(b_gmm,X,y,W0,Z,n,k)

g = 1/n.*Z'*(y-X*b_gmm);
Q = g'*W0*g;

end