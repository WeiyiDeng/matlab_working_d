function fire07try(xmax,ymax,total_percent,spread_speed,move_speed,lambda)

%clear all
global x y r i overlap                           % 设置全局变量，这样在function里面也可以调用这些变量

rand('twister', sum(100*clock))  
axis([0,xmax,0,ymax])
daspect([1,1,1]) 
k=1;d=1;df=0;r(k)=0;x(k)=0;y(k)=0;i=0;
interarrival_day = round(exprnd(lambda))              % 第一个火头的发生时间服从指数分布

for d =1:365         
    if d == df + interarrival_day                
      %if d == df+10
           df=d;
           x_new = unifrnd(0,xmax,1)
           y_new = unifrnd(0,ymax,1)
           interarrival_day = round(exprnd(lambda))    % 下一个火头的发生时间服从指数分布,见fire04try2里的注释
           if i > 1
               judge_overlap(x_new,y_new)       %调用 function judge_overlap 判断新火头的坐标是否与任意旧火范围重合
               if any(overlap)==0               %当与所有旧火都没有重合时，产生新火,见fire04try2里的注释            
                   i=i+1;
                   x(i)=x_new                                
                   y(i)=y_new
                   r(i)=spread_speed
                   %draw_fire(x(i),y(i),r(i))
                   rectangle('Position',[x(i)-r(i),y(i)-r(i),2*r(i),2*r(i)],'Curvature',[1,1],...
     'FaceColor','r') 
               else
               end
           else                                 % 对第一个火头(i=1)无需判断是否重合
               i=i+1;
               x(i)=unifrnd(0,xmax,1) 
               y(i)=unifrnd(0,ymax,1)
               r(i)=spread_speed
           end
      end
      for k=1:i
      x(k)=x(k)+(-1+2*rand)*move_speed 
      y(k)=y(k)+(-1+2*rand)*move_speed
      r(k)=r(k)+spread_speed
      %draw_fire(x(k),y(k),r(k))

   rectangle('Position',[x(k)-r(k),y(k)-r(k),2*r(k),2*r(k)],'Curvature',[1,1],...
     'FaceColor','r')            
      end
      
       pause(0.05)
       
x_space = 0.5:1:(xmax-0.5);                           % 在长度为130km的x轴上，写出每km的中点
y_space = 0.5:1:(ymax-0.5);                           % 在长度为240km的y轴上，写出每km的中点  通过这两步可以确定130*240网格上每一格中心点的坐标
for m = 1:length(x_space)
    for n = 1:length(y_space)
        judge_overlap(m,n)                           %调用function judge_overlap 对每一点判断该点是否在各个火圈的范围内
        if any(overlap)==1                           %只要该点在任意一个火圈的范围内，就算这平方公里着火了
            point_onfire(m,n)=1;
        else
            point_onfire(m,n)=0;                     %否则算没着火
        end
    end
end
area_burned = sum(point_onfire(:) == 1);             %算出总共有多少点（平方公里）着火了
percent_burned = area_burned/(xmax*ymax);              %着火面积所占的比例
if percent_burned > total_percent                             %着火面积大于1%则跳出循环
    break
end
end

day = d                                              %输出在第几天着火面积大于设定比例
      
 
   

  

