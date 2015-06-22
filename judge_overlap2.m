function overlap = judge_overlap2(new_xcoord,new_ycoord, fire_num, x_fire, y_fire, r_fire)        
assert(fire_num >= 1);                              

for j = 1:fire_num
               if (new_xcoord-x_fire(j))^2+(new_ycoord-y_fire(j))^2<=r_fire(j)^2       
                   overlap(j)=1                               
               else
                   overlap(j)=0                              
               end
end

end