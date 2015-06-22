function overlap = judge_overlap(x_coord,y_coord)          % 和fire04try里的内容一样，只是写成了function
global x y r i overlap                                     % 设置全局变量
for j = 1:i
               if (x_coord-x(j))^2+(y_coord-y(j))^2<=r(j)^2       
                   overlap(j)=1                               
               else
                   overlap(j)=0                              
               end
end