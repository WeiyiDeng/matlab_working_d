rand('twister', sum(100*clock))  %初始化随机数状态
axis([0,130,0,240])
day=1
for day=1:360
    while 1  %1为真
       if rem(day,8)==0 %天数除以8余数为0
           break  %跳出while循环
       end
       day=day+1
    end
    x=unifrnd(0,130,1) %在0-130产生一个随机数赋给x
    y=unifrnd(0,240,1)
    R=0.4*day
    h=line(x,y,'color','red','marker','.','markersize',R)
    pause(0.05)
end
  

