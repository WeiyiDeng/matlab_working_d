function fire05(xmax,ymax)

rand('twister', sum(100*clock))  %��ʼ�������״̬
axis([0,xmax,0,ymax])
daspect([1,1,1]) 
k=1;d=1;df=0;r(k)=0;x(k)=0;y(k)=0;i=1;flag_new_fire=0;

for d =1:200
   
     if d == df+round(exprnd(8))
     %if d == df+10
          flag_new_fire=0;  %�»������־��Ĭ��0--�»������1--��������
          x_new=unifrnd(0,xmax,1); %��0-130����һ�����������x
          y_new=unifrnd(0,ymax,1); %��0-240����һ�����������y
         %����ǵ�һ����ͷ��i=1,�����жϣ������ж��Ƿ������л𳡷�Χ�ڣ�
         if i==1 
           df=d;
           x(1)=x_new;  %��һ����ͷ�����ġ��뾶
           y(1)=y_new;
           r(1)=0.4;
           rectangle('Position',[x(1)-r(1),y(1)-r(1),2*r(1),2*r(1)],'Curvature',[1,1],...
          'FaceColor','r')       
         elseif  i>1 
            for j=1:i     %�����л�ͷ����Ƚϣ����»��������κ�һ��j���л�ͷ���ľ���С����뾶rj,���»������־��Ϊ1�����������»�
                if (x_new-x(j))^2+(y_new-y(j))^2<=r(j)^2
                  flag_new_fire=1;
                  % return   %��Ϊ���л�ͷ�������ޣ���������ѭ��Ҳ�����˷�̫��ʱ��
                  % break
                end
            end
         end
        if flag_new_fire==0 
           df=d;i=i+1;
           x(i)=x_new; %��0-130����һ�������������i����ͷ��
           y(i)=y_new;
           r(i)=0.4;
           rectangle('Position',[x(i)-r(i),y(i)-r(i),2*r(i),2*r(i)],'Curvature',[1,1],...
          'FaceColor','r')       
         end
      else %������if d == df+round(exprnd(8))�����������»�������л�ͷ�����ƶ���������
      for k=1:i %iΪ���л�ͷ���� 
      r(k)=r(k)+0.4;        
      x(k)=x(k)-0.6+1.2*rand; %x������ƶ� sqrt(0.85^2/2)=0.6
      y(k)=y(k)-0.6+1.2*rand; %y������ƶ�
      %x0=x(k);y0=y(k); 
      %x(k)=x(k)-0.85+2*0.85*rand;  %x������ƶ�һ������-0.85  +0.85��
      %y(k)=sqrt(0.85*0.85-(x(k)-x0)*(x(k)-x0))+y0;%ͨ������0.85���y(k)
     %solve('(x(k)-x0)^2+(y(k)-y0)^2=0.85*0.85','y(k)')%�ⷽ����y���ƶ�������֤ÿ���ƶ�0.85km
     
   % h=line(x,y,'color','red','marker','.','markersize',R)
   rectangle('Position',[x(k)-r(k),y(k)-r(k),2*r(k),2*r(k)],'Curvature',[1,1],...
     'FaceColor','r')            % w: ��rectangle��Բ
      end
     end     
       pause(0.05) %�ӳ�һ����ʱ����ڹ۲�ͼ�α仯
 end