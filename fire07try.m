function fire07try(xmax,ymax,total_percent,spread_speed,move_speed,lambda)

%clear all
global x y r i overlap                           % ����ȫ�ֱ�����������function����Ҳ���Ե�����Щ����

rand('twister', sum(100*clock))  
axis([0,xmax,0,ymax])
daspect([1,1,1]) 
k=1;d=1;df=0;r(k)=0;x(k)=0;y(k)=0;i=0;
interarrival_day = round(exprnd(lambda))              % ��һ����ͷ�ķ���ʱ�����ָ���ֲ�

for d =1:365         
    if d == df + interarrival_day                
      %if d == df+10
           df=d;
           x_new = unifrnd(0,xmax,1)
           y_new = unifrnd(0,ymax,1)
           interarrival_day = round(exprnd(lambda))    % ��һ����ͷ�ķ���ʱ�����ָ���ֲ�,��fire04try2���ע��
           if i > 1
               judge_overlap(x_new,y_new)       %���� function judge_overlap �ж��»�ͷ�������Ƿ�������ɻ�Χ�غ�
               if any(overlap)==0               %�������оɻ�û���غ�ʱ�������»�,��fire04try2���ע��            
                   i=i+1;
                   x(i)=x_new                                
                   y(i)=y_new
                   r(i)=spread_speed
                   %draw_fire(x(i),y(i),r(i))
                   rectangle('Position',[x(i)-r(i),y(i)-r(i),2*r(i),2*r(i)],'Curvature',[1,1],...
     'FaceColor','r') 
               else
               end
           else                                 % �Ե�һ����ͷ(i=1)�����ж��Ƿ��غ�
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
       
x_space = 0.5:1:(xmax-0.5);                           % �ڳ���Ϊ130km��x���ϣ�д��ÿkm���е�
y_space = 0.5:1:(ymax-0.5);                           % �ڳ���Ϊ240km��y���ϣ�д��ÿkm���е�  ͨ������������ȷ��130*240������ÿһ�����ĵ������
for m = 1:length(x_space)
    for n = 1:length(y_space)
        judge_overlap(m,n)                           %����function judge_overlap ��ÿһ���жϸõ��Ƿ��ڸ�����Ȧ�ķ�Χ��
        if any(overlap)==1                           %ֻҪ�õ�������һ����Ȧ�ķ�Χ�ڣ�������ƽ�������Ż���
            point_onfire(m,n)=1;
        else
            point_onfire(m,n)=0;                     %������û�Ż�
        end
    end
end
area_burned = sum(point_onfire(:) == 1);             %����ܹ��ж��ٵ㣨ƽ������Ż���
percent_burned = area_burned/(xmax*ymax);              %�Ż������ռ�ı���
if percent_burned > total_percent                             %�Ż��������1%������ѭ��
    break
end
end

day = d                                              %����ڵڼ����Ż���������趨����
      
 
   

  

