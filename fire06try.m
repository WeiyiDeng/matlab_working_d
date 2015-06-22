clear all

rand('twister', sum(100*clock))  %初始化随机数状态
axis([0,130,0,240])
daspect([1,1,1]) 
k=1;d=1;df=0;r(k)=0;x(k)=0;y(k)=0;i=0;
for d =1:200   
      
    %if d == df+round(exprnd(8))
      if d == df+10
           df=d;
           x_new = unifrnd(0,130,1)
           y_new = unifrnd(0,240,1)
           if i > 1
               for j = 1:i
               if (x_new-x(j))^2+(y_new-y(j))^2<=r(j)^2       %判断新的着火点是否在第j个火圈的范围内
                   overlap(j)=1                               %把新的着火点座标与每个火圈的比较结果存在向量overlap内，向量overlap由i个0-1变量组成
               else
                   overlap(j)=0                               %如果新的着火点在第j个火圈的范围内，给向量overlap的第j个元素赋值1，否则为0
               end
           end
               if any(overlap)==0                             %如果向量overlap的全部组成元素都是0，即新的着火点与现有的任何火圈不重合，产生新火苗
                   i=i+1;
                   x(i)=x_new                                 %新起火点的座标是之前产生的随机数
                   y(i)=y_new
                   r(i)=0.4
                   rectangle('Position',[x(i)-r(i),y(i)-r(i),2*r(i),2*r(i)],'Curvature',[1,1],...
     'FaceColor','r') 
               else
               end
           else                      %如果是第一个火头即i=1,不需判断
               i=i+1;
               x(i)=unifrnd(0,130,1) %在0-130产生一个随机数赋给x
               y(i)=unifrnd(0,240,1)
               r(i)=0.4
           end
      end
      for k=1:i
      x(k)=x(k)-1+2*rand %x点随机移动
      y(k)=y(k)-1+2*rand %x点随机移动
      r(k)=r(k)+0.4
   % h=line(x,y,'color','red','marker','.','markersize',R)
   rectangle('Position',[x(k)-r(k),y(k)-r(k),2*r(k),2*r(k)],'Curvature',[1,1],...
     'FaceColor','r')            % w: 用rectangle画圆
      end
      
       pause(0.05)
end
      
      
 
   

  

