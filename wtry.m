clear all
global x y r i overlap                           % 设置全局变量，这样在function里面也可以调用这些变量

rand('twister', sum(100*clock))  
axis([0,130,0,240])
daspect([1,1,1]) 
k=1;d=1;df=0;r(k)=0;x(k)=0;y(k)=0;i=0;
interarrival_day = round(exprnd(8))              % 第一个火头的发生时间服从指数分布

for d =1:200         
    if d == df + interarrival_day                
      %if d == df+10
           df=d;
           x_new = unifrnd(0,130,1)
           y_new = unifrnd(0,240,1)
           interarrival_day = round(exprnd(8))    % 下一个火头的发生时间服从指数分布,见fire04try2里的注释
           if i > 1
               judge_overlap(x_new,y_new)       %调用 function judge_overlap 判断新火头的坐标是否与任意旧火范围重合
               if any(overlap)==0               %当与所有旧火都没有重合时，产生新火,见fire04try2里的注释            
                   i=i+1;
                   x(i)=x_new                                
                   y(i)=y_new
                   r(i)=0.4
                   rectangle('Position',[x(i)-r(i),y(i)-r(i),2*r(i),2*r(i)],'Curvature',[1,1],...
     'FaceColor','r') 
               else
               end
           else                                 % 对第一个火头(i=1)无需判断是否重合
               i=i+1;
               x(i)=unifrnd(0,130,1) 
               y(i)=unifrnd(0,240,1)
               r(i)=0.4
           end
      end
      for k=1:i
      x(k)=x(k)-1+2*rand 
      y(k)=y(k)-1+2*rand
      r(k)=r(k)+0.4

   rectangle('Position',[x(k)-r(k),y(k)-r(k),2*r(k),2*r(k)],'Curvature',[1,1],...
     'FaceColor','r')            
      end
      
       pause(0.05)
       

end

day = d                                              %输出在第几天着火面积大于设定比例
      
 
   

  

