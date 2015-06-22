rand('twister', sum(100*clock))  %初始化随机数状态
axis([0,130,0,240])

for i = 1:32
    if rem(i,8)==0
        x(i)=unifrnd(0,130,1) %在0-130产生一个随机数赋给x
        y(i)=unifrnd(0,240,1)
        num(i) = i
    else        
    end
    for j = 1:i
        if x(j) > 0
            R(j) = 0.4*(i-num(j))
        else
            R(j)=0
        end
    end
    if R(i) > 0
    rectangle('Position',[x(i)-R(i),y(i)-R(i),2*R(i),2*R(i)],'Curvature',[1,1],...
     'FaceColor','r')            % w: 用rectangle画圆
    pause(0.05)
    else
    end
end
    