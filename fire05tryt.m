clear all
%global x y r i overlap                  No need for global variable now

rand('twister', sum(100*clock))  
axis([0,130,0,240])
daspect([1,1,1]) 
k=1;d=1;df=0;r(k)=0;x(k)=0;y(k)=0;i=0;
interarrival_day = round(exprnd(8))              

for d =1:365
    cla                                 % clean up the track of move
    if d == df + interarrival_day                
      %if d == df+10
           df=d;
           x_new = unifrnd(0,130,1)
           y_new = unifrnd(0,240,1)
           interarrival_day = round(exprnd(8))    
           if i > 1
               retval = judge_overlap2(x_new,y_new, i, x, y, r)        % store the result of judge_overlap in retval
               if any(retval)==0                      
                   i=i+1;
                   x(i)=x_new                                
                   y(i)=y_new
                   r(i)=0.4
                   rectangle('Position',[x(i)-r(i),y(i)-r(i),2*r(i),2*r(i)],'Curvature',[1,1],...
     'FaceColor','r') 
               else
               end
           else                                 
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
       
x_space = 0.5:1:(130-0.5);                         
y_space = 0.5:1:(240-0.5);
%point_onfire=zeros(300,500,'single');
if i >= 1
for m = 1:length(x_space)
    for n = 1:length(y_space)
        retval = judge_overlap2(m,n,i, x, y, r)                       
        if any(retval)==1                           
            point_onfire(m,n)=1;
        else
            point_onfire(m,n)=0;                   
        end
    end
end
area_burned = sum(point_onfire(:) == 1);            
percent_burned = area_burned/(130*240);             
if percent_burned > 0.1                            
    break
end
else
end
end

day = d                                           
      
 
   

  

