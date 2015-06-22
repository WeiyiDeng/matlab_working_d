for i = 1:360
    if rem(i,8)==0
        x(i)=unifrnd(0,130,1) %在0-130产生一个随机数赋给x
        y(i)=unifrnd(0,240,1)
        firetime(i) = i
    else        
    end
    %for j = 1:i
            %R(j) = 0.4*(i-num(j))
        %end
end

for j = 1:360
    for i = 1:j
        if x(i) > 0 
        R(i) = 0.4*(j - firetime(i))
        else
        R(i) = 0
        end
    end
end
