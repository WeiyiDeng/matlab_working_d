clear all
%global x y r i overlap                  No need for global variables now

rand('twister', sum(100*clock))  
axis([0,130,0,240])
daspect([1,1,1]) 
%k=1;r(k)=0;x(k)=0;y(k)=0;
day=1;day_formerfire=0;
fire_index=0;
interarrival_day = round(exprnd(8))              

for day =1:365
    cla                                 % clean up the track of move
     if day == day_formerfire + interarrival_day                 
      %if d == df+10
           day_formerfire = day;
           x_newfire = unifrnd(0,130,1)
           y_newfire = unifrnd(0,240,1)
           interarrival_day = round(exprnd(8))    
           if fire_index > 1
               retval = check_burned(x_newfire,y_newfire, fire_index, xcoord_fire, ycoord_fire, radius_fire)        % store the result of judge_overlap in retval
               if any(retval)==0                      
                   fire_index=fire_index+1;
                   xcoord_fire(fire_index)=x_newfire                                
                   ycoord_fire(fire_index)=y_newfire
                   radius_fire(fire_index)=0.4
                  % plot_fire(i, x(i), y(i), r(i))
                 rectangle('Position',[xcoord_fire(fire_index)-radius_fire(fire_index),ycoord_fire(fire_index)-radius_fire(fire_index),2*radius_fire(fire_index),2*radius_fire(fire_index)],'Curvature',[1,1],...
     'FaceColor','r') 
               else
               end
           else                                 
               fire_index=fire_index+1;
               xcoord_fire(fire_index)=unifrnd(0,130,1) 
               ycoord_fire(fire_index)=unifrnd(0,240,1)
               radius_fire(fire_index)=0.4
           end
      end
      for k=1:fire_index
      xcoord_fire(k)=xcoord_fire(k)-1+2*rand 
      ycoord_fire(k)=ycoord_fire(k)-1+2*rand
      radius_fire(k)=radius_fire(k)+0.4
     % plot_fire(k, x(k), y(k), r(k))
   rectangle('Position',[xcoord_fire(k)-radius_fire(k),ycoord_fire(k)-radius_fire(k),2*radius_fire(k),2*radius_fire(k)],'Curvature',[1,1],...
     'FaceColor','r')            
      end
      
     pause(0.01)
       
x_space = 0.5:1:(130-0.5);                         
y_space = 0.5:1:(240-0.5);
%point_onfire=zeros(300,500,'single');
if fire_index >= 1
for m = 1:length(x_space)                     % compare each point on the grid with the center of each fire circle
    for n = 1:length(y_space)
        retval = check_burned(m,n,fire_index, xcoord_fire, ycoord_fire, radius_fire)                       
        if any(retval)==1                           
            point_burned(m,n)=1;
        else
            point_burned(m,n)=0;                   
        end
    end
end
area_burned = sum(point_burned(:) == 1);            
percent_burned = area_burned/(130*240);             
if percent_burned > 0.1                            
    break
end
else
end
end

day                                           
      
 
   

  

