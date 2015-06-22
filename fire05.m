function fire05(xmax,ymax)

rand('twister', sum(100*clock))  %初始化随机数状态
axis([0,xmax,0,ymax])
daspect([1,1,1]) 
k=1;d=1;df=0;r(k)=0;x(k)=0;y(k)=0;i=1;flag_new_fire=0;

for d =1:200
   
     if d == df+round(exprnd(8))
     %if d == df+10
          flag_new_fire=0;  %新火产生标志，默认0--新火产生，1--不产生；
          x_new=unifrnd(0,xmax,1); %在0-130产生一个随机数赋给x
          y_new=unifrnd(0,ymax,1); %在0-240产生一个随机数赋给y
         %如果是第一个火头即i=1,则不需判断；否则判断是否在已有火场范围内；
         if i==1 
           df=d;
           x(1)=x_new;  %第一个火头的中心、半径
           y(1)=y_new;
           r(1)=0.4;
           rectangle('Position',[x(1)-r(1),y(1)-r(1),2*r(1),2*r(1)],'Curvature',[1,1],...
          'FaceColor','r')       
         elseif  i>1 
            for j=1:i     %与已有火头逐个比较，若新火中心与任何一个j已有火头中心距离小于其半径rj,则将新火产生标志改为1（即不产生新火）
                if (x_new-x(j))^2+(y_new-y(j))^2<=r(j)^2
                  flag_new_fire=1;
                  % return   %因为已有火头总数有限，不用跳出循环也不会浪费太多时间
                  % break
                end
            end
         end
        if flag_new_fire==0 
           df=d;i=i+1;
           x(i)=x_new; %在0-130产生一个随机数赋给第i个火头；
           y(i)=y_new;
           r(i)=0.4;
           rectangle('Position',[x(i)-r(i),y(i)-r(i),2*r(i),2*r(i)],'Curvature',[1,1],...
          'FaceColor','r')       
         end
      else %不满足if d == df+round(exprnd(8))，即不产生新火，则对已有火头进行移动、扩大处理
      for k=1:i %i为已有火头总数 
      r(k)=r(k)+0.4;        
      x(k)=x(k)-0.6+1.2*rand; %x点随机移动 sqrt(0.85^2/2)=0.6
      y(k)=y(k)-0.6+1.2*rand; %y点随机移动
      %x0=x(k);y0=y(k); 
      %x(k)=x(k)-0.85+2*0.85*rand;  %x点随机移动一个量（-0.85  +0.85）
      %y(k)=sqrt(0.85*0.85-(x(k)-x0)*(x(k)-x0))+y0;%通过距离0.85解得y(k)
     %solve('(x(k)-x0)^2+(y(k)-y0)^2=0.85*0.85','y(k)')%解方程求y点移动量，保证每天移动0.85km
     
   % h=line(x,y,'color','red','marker','.','markersize',R)
   rectangle('Position',[x(k)-r(k),y(k)-r(k),2*r(k),2*r(k)],'Curvature',[1,1],...
     'FaceColor','r')            % w: 用rectangle画圆
      end
     end     
       pause(0.05) %延迟一定的时间便于观察图形变化
 end