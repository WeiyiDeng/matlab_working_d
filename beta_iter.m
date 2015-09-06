function  [stop,options,optchanged] = beta_iter(options,optimvalues,flag)
% function state = beta_iter(b,options,state)
global b_best ll_best

b = optimvalues.bestx;
ll = optimvalues.bestfval;
% if isequal(flag,'iter')
%     b_hist = [b_hist; b];
%     ll_hist = [ll_hist; ll];
% end
if ll < ll_best
    ll_best = ll;
    b_best = b;
    save('b_best.mat','b_best','-append')
    save('ll_best.mat','ll_best','-append')
end

stop = false;
optchanged = false;
% switch flag
%    case 'init'
%         disp('Initializing output function');
%     case 'iter'
%         disp('Iterating ...')
%     case 'done'
%         disp('Performing final task');

end

% check out the template output function file
% edit saoutputfcntemplate