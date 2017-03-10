%% generate data
rng(10086)

% sth = ceil(abs(randn(1,20)*10))
myarray = ceil(abs(trnd(1,1,20)*5-3))

plot(myarray,'o')

%%
k=0;
smooth_array = zeros(size(myarray));

for t=1:length(myarray)
    if myarray(t)>10
        break
    end
end
if t > 1
    smooth_array(1:(t-1)) = mean(myarray(1:t));
else
end
begin_ind = t

for i=begin_ind:length(myarray)
    if myarray(i)>10
        k = 0;
        smooth_array(i) = myarray(i);
    else k = k+1;
        right_end_add = 0;
        for j = 1:k
            if i+j > length(myarray)
                break
            elseif myarray(i+j)>10
                right_end_add = j;
                break
            end
        end
        if right_end_add == 0 && (i+j < length(myarray))
            smooth_array(i) = sum(myarray((i-k):(i+k)))/(2*k+1);
        elseif (i+j < length(myarray))
            smooth_array(i) = sum(myarray((i-j):(i+j)))/(2*j+1);
        else smooth_array(i) = 999;           % for indices exceeding length(myarray), to change the value manually later
        end
    end
end

% Index exceeds matrix dimensions            to be fixed ?

% % test break
% for i = 1:10
%     for j = 1:10
%         my_times = i*j
%         if my_times >50
%             break                 % only breaks the first loop
%         end
%     end
% end