function [z, num_itera] = itera(z,c,n)
z = z^2 + c;
num_itera = n;
if abs(z) > 2
else n = n + 1;
    if n > 100
        num_itera = 101;
    else
    [z, num_itera] = itera(z,c,n);
    end
end;

