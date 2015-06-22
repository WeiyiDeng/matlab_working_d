function plot_fire(fire_index, x_fire, y_fire ,r_fire)

rectangle('Position',[x_fire(fire_index)-r_fire(fire_index),y_fire(fire_index)-r_fire(fire_index),2*r_fire(fire_index),2*r_fire(fire_index)],'Curvature',[1,1],...
     'FaceColor','r')
end