function overlap = judge_overlap(x_coord,y_coord)          % ��fire04try�������һ����ֻ��д����function
global x y r i overlap                                     % ����ȫ�ֱ���
for j = 1:i
               if (x_coord-x(j))^2+(y_coord-y(j))^2<=r(j)^2       
                   overlap(j)=1                               
               else
                   overlap(j)=0                              
               end
end