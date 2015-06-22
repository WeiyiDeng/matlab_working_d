clear all

rand('twister', sum(100*clock))  
axis([0,130,0,240])
daspect([1,1,1]) 
k=1;d=1;df=0;r(k)=0;x(k)=0;y(k)=0;i=0;
interarrival_day = round(exprnd(8))           

for d =1:200   
    cla                                           % 清除之前移动痕迹
     if d == df + interarrival_day
     %if d == df+10
           df=d;
           x_new = unifrnd(0,130,1)
           y_new = unifrnd(0,240,1)
           interarrival_day = round(exprnd(8))               
           if i > 1
               for j = 1:i
               if (x_new-x(j))^2+(y_new-y(j))^2<=r(j)^2       
                   overlap(j)=1                              
               else
                   overlap(j)=0                              
               end
           end
               if any(overlap)==0                           
                   i=i+1;
                   x(i)=x_new                               
                   y(i)=y_new
                   r(i)=0.4
                  
                   rectangle('Position',[x(i)-r(i),y(i)-r(i),2*r(i),2*r(i)],'Curvature',[1,1],...
     'FaceColor','r') 
               else
               end
           else                      
               i=i+1;
               x(i)=unifrnd(0,130,1) 
               y(i)=unifrnd(0,240,1)
               r(i)=0.4
           end
      end
      for k=1:i
      x(k)=x(k)+(-1+2*rand)*5                      %随机移动速度乘5
      y(k)=y(k)+(-1+2*rand)*5                      %随机移动速度乘5
      r(k)=r(k)+0.4
   % h=line(x,y,'color','red','marker','.','markersize',R)
   rectangle('Position',[x(k)-r(k),y(k)-r(k),2*r(k),2*r(k)],'Curvature',[1,1],...
     'FaceColor','r')            
      end
      
       pause(0.05)
end
      
      
 
   

  

