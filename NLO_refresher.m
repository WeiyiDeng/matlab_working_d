%% w: my approach  (better approach !!)
% x1 = linspace(-1,3,100);
% x2 = linspace(-1,3,100);
% f = zeros(length(x1),length(x2));
% for i = 1:length(x1)
% for j = 1:length(x2)
% f(i,j) = 100*(x2(j)-x1(i)^2)^2+(1-x1(i))^2;
% end
% end
% surfc(x1,x2,f)

x1 = linspace(-2,2,100);
x2 = linspace(-2,2,100);
z = zeros(length(x1),length(x2));
f = @(x,y) 100*(x-y^2)^2+(1-x)^2;             % Rosenbrock(banana) function
for i = 1:length(x1)
    for j = 1:length(x2)
        z(i,j) = f(x1(i),x2(j));
    end
end
surfc(x1,x2,z)                                % input + x1 & x2 to give the correct coordinates!!
xlabel('x');
ylabel('y');
zlabel('f(x,y) = Rosenbrock function');
hold on
plot3(1, 1, f(1,1),'.r','MarkerSize',30);       % mark the minimun point
hold off

% figure(gcf)                                   % return current plot        

%% copy from course material
x = linspace(-1,1);
y = linspace(-1,1);
[X, Y] = meshgrid(x,y);                       % w: no need to use meshgrid

f = @(x,y) x*y; %define function f
Z = zeros(size(X)); %initialize matrix Z
for i = 1:size(X,1) %use a double loop to calculate all values
for j = 1:size(X,2)
Z(i,j) = f(X(i,j),Y(i,j));
end
end
% f = @(X,Y) X.*Y; %use .* to define elementwise multiplication
% Z = f(X,Y); %determine Z in one go

surfc(X,Y,Z);
xlabel('x');
ylabel('y');
zlabel('f(x,y) = xy');

hold on %enable working in the same figure
plot3(0.5, 0.5, f(0.5,0.5),'.r','MarkerSize',30);
hold off

figure(gcf)
