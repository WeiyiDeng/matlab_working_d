rand('twister', sum(100*clock))  %��ʼ�������״̬
axis([0,130,0,240])
day=1
for day=1:360
    while 1  %1Ϊ��
       if rem(day,8)==0 %��������8����Ϊ0
           break  %����whileѭ��
       end
       day=day+1
    end
    x=unifrnd(0,130,1) %��0-130����һ�����������x
    y=unifrnd(0,240,1)
    R=0.4*day
    h=line(x,y,'color','red','marker','.','markersize',R)
    pause(0.05)
end
  

