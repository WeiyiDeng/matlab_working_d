clear all

rand('twister', sum(100*clock))  %��ʼ�������״̬
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
               if (x_new-x(j))^2+(y_new-y(j))^2<=r(j)^2       %�ж��µ��Ż���Ƿ��ڵ�j����Ȧ�ķ�Χ��
                   overlap(j)=1                               %���µ��Ż��������ÿ����Ȧ�ıȽϽ����������overlap�ڣ�����overlap��i��0-1�������
               else
                   overlap(j)=0                               %����µ��Ż���ڵ�j����Ȧ�ķ�Χ�ڣ�������overlap�ĵ�j��Ԫ�ظ�ֵ1������Ϊ0
               end
           end
               if any(overlap)==0                             %�������overlap��ȫ�����Ԫ�ض���0�����µ��Ż�������е��κλ�Ȧ���غϣ������»���
                   i=i+1;
                   x(i)=x_new                                 %�������������֮ǰ�����������
                   y(i)=y_new
                   r(i)=0.4
                   rectangle('Position',[x(i)-r(i),y(i)-r(i),2*r(i),2*r(i)],'Curvature',[1,1],...
     'FaceColor','r') 
               else
               end
           else                      %����ǵ�һ����ͷ��i=1,�����ж�
               i=i+1;
               x(i)=unifrnd(0,130,1) %��0-130����һ�����������x
               y(i)=unifrnd(0,240,1)
               r(i)=0.4
           end
      end
      for k=1:i
      x(k)=x(k)-1+2*rand %x������ƶ�
      y(k)=y(k)-1+2*rand %x������ƶ�
      r(k)=r(k)+0.4
   % h=line(x,y,'color','red','marker','.','markersize',R)
   rectangle('Position',[x(k)-r(k),y(k)-r(k),2*r(k),2*r(k)],'Curvature',[1,1],...
     'FaceColor','r')            % w: ��rectangle��Բ
      end
      
       pause(0.05)
end
      
      
 
   

  

