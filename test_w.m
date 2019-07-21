x = -1:0.1:1
y = -1:0.1:1
x_b = -1
y_b = -1
x_minus_y_b = 0
x_times_y_b = 3
x_quad_b = 0
y_quad_b = 0
c = 0

f = zeros(length(x),length(y));
for i = 1:length(x)
for j = 1:length(y)
% f(i,j) = sqrt(x_b*x(i)+ y_b*y(j)+x_minus_y_b*(x(i)-y(j))+x_times_y_b*x(i)*y(j)...
%     +x_squared_b*(x(i)^2)+y_squared_b*(y(j)^2))+c;
f(i,j) = x_b*x(i)+ y_b*y(j)+x_minus_y_b*(x(i)-y(j))+x_times_y_b*x(i)*y(j)...
    +x_quad_b*(x(i)^2)+y_quad_b*(y(j)^2)+c;
end
end
surf(f)


% x = -10:1:10
% y = -10:1:10
% x_quad_b = 1
% y_quad_b = 1
% c= 0
% f = zeros(length(x),length(y));
% for i = 1:length(x)
% for j = 1:length(y)
% % f(i,j) = sqrt(x_quad_b*(x(i)^2)+y_quad_b*(y(j)^2))+c;
% f(i,j) = sqrt(x(i)^2+y(j)^2)+1;
% end
% end
% surf(f)

