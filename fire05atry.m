clear all
global x y r i overlap                           % ����ȫ�ֱ�����������function����Ҳ���Ե�����Щ����

rand('twister', sum(100*clock))  
axis([0,130,0,240])
daspect([1,1,1]) 
k=1;d=0;df=0;r(k)=0;x(k)=0;y(k)=0;i=0; percent_burned=0;
interarrival_day = round(exprnd(8));              % ��һ����ͷ�ķ���ʱ�����ָ���ֲ�
length_xpoint = length(0.5:1:(130-0.5));
length_ypoint = length(0.5:1:(240-0.5));

while percent_burned < 0.1    %���������<10%,���̼�����
%for d =1:365
  d=d+1;
  if d>=365  %�����������һ�ֹ꣬ͣ����ʾ 
      warndlg('ȫ���޻�','�ﵽ365��')
      break
  end      
  %if d == df + interarrival_day                
  if d == df+8
           df=d;
           x_new = unifrnd(0,130,1);
           y_new = unifrnd(0,240,1);
           interarrival_day = round(exprnd(8));    % ��һ����ͷ�ķ���ʱ�����ָ���ֲ�,��fire04try2���ע��
           if i > 1
               judge_overlap(x_new,y_new)       %���� function judge_overlap �ж��»�ͷ�������Ƿ�������ɻ�Χ�غ�
               if any(overlap)==0               %�������оɻ�û���غ�ʱ�������»�,��fire04try2���ע��            
                   i=i+1;
                   x(i)=x_new;                                
                   y(i)=y_new;
                   r(i)=0.4;
                   rectangle('Position',[x(i)-r(i),y(i)-r(i),2*r(i),2*r(i)],'Curvature',[1,1],...
     'FaceColor','r') 
               else
               end
  else                                 % �Ե�һ����ͷ(i=1)�����ж��Ƿ��غ�
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
      
    pause(0.005)                                %��Ҫ��ͣ
       
   %x_space = 0.5:1:(130-0.5);                      % �ڳ���Ϊ130km��x���ϣ�д��ÿkm���е�
   %y_space = 0.5:1:(240-0.5);                      % �ڳ���Ϊ240km��y���ϣ�д��ÿkm���е�  ͨ������������ȷ��130*240������ÿһ�����ĵ������
   point_onfire=zeros(300,500,'single');            % Ԥ��0������ٶ�
if d > 45
  for m = 1:length_xpoint
      for n = 1:length_xpoint
        judge_overlap(m,n)                         %����function judge_overlap ��ÿһ���жϸõ��Ƿ��ڸ�����Ȧ�ķ�Χ��
        if any(overlap)==1                         %ֻҪ�õ�������һ����Ȧ�ķ�Χ�ڣ�������ƽ�������Ż���
            point_onfire(m,n)=1;
        else
            point_onfire(m,n)=0;                     %������û�Ż�
        end
      end
  end
else
end
  area_burned = sum(point_onfire(:) == 1);             %����ܹ��ж��ٵ㣨ƽ������Ż���
  percent_burned = area_burned/(130*240);              %�Ż������ռ�ı���
  if percent_burned >= 0.1                             %�Ż��������10%������ѭ��
    warndlg('�����������10%','����')
    break
  end
  
end

day = d;                                              %����ڵڼ����Ż���������趨����
      
 
   

  

