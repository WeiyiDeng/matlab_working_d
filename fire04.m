clear all

rand('twister', sum(100*clock))  %��ʼ�������״̬
axis([0,130,0,240])
daspect([1,1,1]) 
k=1;d=1;df=0;r(k)=0;x(k)=0;y(k)=0;i=0;
for d =1:200
   
      %if d == df+round(exprnd(8))
      if d == df+10
           df=d;i=i+1;
    x(i)=unifrnd(0,130,1) %��0-130����һ�����������x
    y(i)=unifrnd(0,240,1)
    r(i)=0.4
    rectangle('Position',[x(i)-r(i),y(i)-r(i),2*r(i),2*r(i)],'Curvature',[1,1],...
     'FaceColor','r') 
       else
      for k=1:i     
      x(k)=x(k)-1+2*rand %x������ƶ�
      y(k)=y(k)-1+2*rand %x������ƶ�
      r(k)=r(k)+0.4
   % h=line(x,y,'color','red','marker','.','markersize',R)
   rectangle('Position',[x(k)-r(k),y(k)-r(k),2*r(k),2*r(k)],'Curvature',[1,1],...
     'FaceColor','r')            % w: ��rectangle��Բ
      end
      end
       pause(0.05)
   end
      
      
 
   

  

