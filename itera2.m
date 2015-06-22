function [z, num_itera] = itera2(z,c,n)
z = z^2 + c;
%global num_itera
num_itera = n;
if abs(z) > 2;
    n
else n = n + 1;
    if n > 100
        num_itera = 101;
    else
    [z, num_itera] = itera2(z,c,n);
    end
end;

