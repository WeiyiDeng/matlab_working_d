clear all
global x y r i overlap                           % 设置全局变量，这样在function里面也可以调用这些变量

rand('twister', sum(100*clock))  
axis([0,130,0,240])
daspect([1,1,1]) 
k=1;d=0;df=0;r(k)=0;x(k)=0;y(k)=0;i=0; percent_burned=0;
interarrival_day = round(exprnd(8));              % 第一个火头的发生时间服从指数分布
length_xpoint = length(0.5:1:(130-0.5));
length_ypoint = length(0.5:1:(240-0.5));

while percent_burned < 0.1    %若过火面积<10%,进程继续；
%for d =1:365
  d=d+1;
  if d>=365  %如果天数多于一年，停止且提示 
      warndlg('全年无火警','达到365天')
      break
  end      
  %if d == df + interarrival_day                
  if d == df+8
           df=d;
           x_new = unifrnd(0,130,1);
           y_new = unifrnd(0,240,1);
           interarrival_day = round(exprnd(8));    % 下一个火头的发生时间服从指数分布,见fire04try2里的注释
           if i > 1
               judge_overlap(x_new,y_new)       %调用 function judge_overlap 判断新火头的坐标是否与任意旧火范围重合
               if any(overlap)==0               %当与所有旧火都没有重合时，产生新火,见fire04try2里的注释            
                   i=i+1;
                   x(i)=x_new;                                
                   y(i)=y_new;
                   r(i)=0.4;
                   rectangle('Position',[x(i)-r(i),y(i)-r(i),2*r(i),2*r(i)],'Curvature',[1,1],...
     'FaceColor','r') 
               else
               end
  else                                 % 对第一个火头(i=1)无需判断是否重合
               i=i+1;
               x(i)=unifrnd(0,130,1); 
               y(i)=unifrnd(0,240,1);
               r(i)=0.4;
           end
      end
      for k=1:i
      x(k)=x(k)-0.6+1.2*rand; 
      y(k)=y(k)-0.6+1.2*rand;
      r(k)=r(k)+0.4;

   rectangle('Position',[x(k)-r(k),y(k)-r(k),2*r(k),2*r(k)],'Curvature',[1,1],...
     'FaceColor','r')            
      end
      
    pause(0.005)                                %不要暂停
       
   %x_space = 0.5:1:(130-0.5);                      % 在长度为130km的x轴上，写出每km的中点
   %y_space = 0.5:1:(240-0.5);                      % 在长度为240km的y轴上，写出每km的中点  通过这两步可以确定130*240网格上每一格中心点的坐标
   point_onfire=zeros(300,500,'single');            % 预置0，提高速度
if d > 45
  for m = 1:length_xpoint
      for n = 1:length_xpoint
        judge_overlap(m,n)                         %调用function judge_overlap 对每一点判断该点是否在各个火圈的范围内
        if any(overlap)==1                         %只要该点在任意一个火圈的范围内，就算这平方公里着火了
            point_onfire(m,n)=1;
        else
            point_onfire(m,n)=0;                     %否则算没着火
        end
      end
  end
else
end
  area_burned = sum(point_onfire(:) == 1);             %算出总共有多少点（平方公里）着火了
  percent_burned = area_burned/(130*240);              %着火面积所占的比例
  if percent_burned >= 0.1                             %着火面积大于10%则跳出循环
    warndlg('过火面积超过10%','报警')
    break
  end
  
end

day = d;                                              %输出在第几天着火面积大于设定比例
      
 
   

  

