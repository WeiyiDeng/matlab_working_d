rand('twister', sum(100*clock))  %初始化随机数状态
axis([0,130,0,240])
daspect([1,1,1])                 % w: 横纵座标相同比例

x=unifrnd(0,130,1) %在0-130产生一个随机数赋给x
y=unifrnd(0,240,1)

for day=1:360
R=0.4*day
% h=line(x,y,'color','red','marker','.','markersize',R)
rectangle('Position',[x-R,y-R,2*R,2*R],'Curvature',[1,1],...
     'FaceColor','r')            % w: 用rectangle画圆

pause(0.05)
end

 