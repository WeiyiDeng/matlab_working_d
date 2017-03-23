% baseline probability smoothing transformation with ad hoc method 
% (varying moving average window depending on # of observations in the time interval)
function smooth_array = baseline_prob_smooth_func(myarray, obs_cut)

k=0;
smooth_array = zeros(size(myarray));

for t=1:length(myarray)
    if myarray(t)>obs_cut
        break
    end
end
if t > 1
    smooth_array(1:(t-1)) = mean(myarray(1:t));
else
end
begin_ind = t;

for i=begin_ind:length(myarray)
    if myarray(i)>obs_cut            % number of observations used as cut-off point; not to do smoothing if above this point
        k = 0;
        smooth_array(i) = myarray(i);
    else k = k+1;
        right_end_add = 0;
        for j = 1:k
            if i+j > length(myarray)
                break
            elseif myarray(i+j)>obs_cut
                right_end_add = j;
                break
            end
        end
        if right_end_add == 0 && (i+j < length(myarray))
            smooth_array(i) = sum(myarray((i-k):(i+k)))/(2*k+1);
        elseif (i+j < length(myarray))
            smooth_array(i) = sum(myarray((i-j):(i+j)))/(2*j+1);
%         else smooth_array(i) = 999;           % for indices exceeding length(myarray), to change the value manually later
        else smooth_array(i) = myarray(i);
        end
    end
end

end