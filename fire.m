rand('twister', sum(100*clock))  %��ʼ�������״̬
axis([0,130,0,240])
daspect([1,1,1])                 % w: ����������ͬ����

x=unifrnd(0,130,1) %��0-130����һ�����������x
y=unifrnd(0,240,1)

for day=1:360
R=0.4*day
% h=line(x,y,'color','red','marker','.','markersize',R)
rectangle('Position',[x-R,y-R,2*R,2*R],'Curvature',[1,1],...
     'FaceColor','r')            % w: ��rectangle��Բ

pause(0.05)
end

 