clear all
global x y r i overlap                           % ����ȫ�ֱ�����������function����Ҳ���Ե�����Щ����

rand('twister', sum(100*clock))  
axis([0,130,0,240])
daspect([1,1,1]) 
k=1;d=1;df=0;r(k)=0;x(k)=0;y(k)=0;i=0;
interarrival_day = round(exprnd(8))              % ��һ����ͷ�ķ���ʱ�����ָ���ֲ�

for d =1:200         
    if d == df + interarrival_day                
      %if d == df+10
           df=d;
           x_new = unifrnd(0,130,1)
           y_new = unifrnd(0,240,1)
           interarrival_day = round(exprnd(8))    % ��һ����ͷ�ķ���ʱ�����ָ���ֲ�,��fire04try2���ע��
           if i > 1
               judge_overlap(x_new,y_new)       %���� function judge_overlap �ж��»�ͷ�������Ƿ�������ɻ�Χ�غ�
               if any(overlap)==0               %�������оɻ�û���غ�ʱ�������»�,��fire04try2���ע��            
                   i=i+1;
                   x(i)=x_new                                
                   y(i)=y_new
                   r(i)=0.4
                   rectangle('Position',[x(i)-r(i),y(i)-r(i),2*r(i),2*r(i)],'Curvature',[1,1],...
     'FaceColor','r') 
               else
               end
           else                                 % �Ե�һ����ͷ(i=1)�����ж��Ƿ��غ�
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

day = d                                              %����ڵڼ����Ż���������趨����
      
 
   

  
