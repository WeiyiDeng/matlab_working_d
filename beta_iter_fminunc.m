% w: wasted, not working

% % function  [stop,options,optchanged] = beta_iter_fminunc(options,optimvalues,flag)
% % function state = beta_iter(b,options,state)
% function stop = beta_iter_fminunc(x,optimvalues,state)
% global b_best ll_best
% 
% stop = false;
% % optchanged = false;
% % switch flag
% %    case 'init'
% %         disp('Initializing output function');
% %     case 'iter'
% %         disp('Iterating ...')
% %     case 'done'
% %         disp('Performing final task');
% 
% switch state
%     case 'init'
%         disp('Initializing output function');
%     case 'iter'
%         b = x;
%         ll = optimvalues.fval;
%         
%         if ll < ll_best
%             ll_best = ll;
%             b_best = b;
%             save('b_fminunc.mat','b_best','-append')
%             save('ll_fminunc.mat','ll_best','-append')
%         end
%     case 'done'
%         disp('Performing final task');
%     otherwise
% end
% 
% end
% 
% % check out the template output function file
% % edit saoutputfcntemplate