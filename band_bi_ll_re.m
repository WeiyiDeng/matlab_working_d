function [LL_sum, gr, H] = band_bi_ll_re(b, member, friendlist_updated, EA_col)
% global friendlist_updated
% global I J choice_dv IVs se T

% const = repmat(b(1),I*T,1);
LL_sum = 0;
mf_rows = find(friendlist_updated(:,1)==member);
for i = 1:size(mf_rows)     % tbchanged
    % load(sprintf('%d_%d',2898,3))
    % friend_files = cell(3,1);
    % for i = 1:3
    %     friend_files{i} = load(sprintf('C:\\Users\\etp21998@eur.nl\\matfolder\\%d_%d.mat',2898,i));
    % end
    % friend_files{1}.mat_r(1,2)
%     load(sprintf('E:\\Wei\\Graduate\\Matlab\\matfolder\\%d_%d',member,mf_rows(i)));
    load(sprintf('C:\\Users\\etp21998@eur.nl\\matfolder\\%d_%d.mat',member,mf_rows(i)));
%     mat_dummies = [mat_r(:,11).*mat_r(:,15) mat_r(:,11).*mat_r(:,16) mat_r(:,12).*mat_r(:,14) mat_r(:,12).*mat_r(:,15) mat_r(:,12).*mat_r(:,16) mat_r(:,13).*mat_r(:,14) mat_r(:,13).*mat_r(:,15) mat_r(:,13).*mat_r(:,16)];
%     mat_r = [mat_r mat_dummies];
    
    IVs = mat_r(:,[6:8 10 EA_col]);
    choice_dv = [mat_r(:,5) 1-mat_r(:,5)];
    clearvars mat_r
    
    const = b(1);
    
    bs = b(2:end)';
    FV = IVs*bs;
    
    % exp_util = exp(const+FV);          % utility of choosing the product
    % prob=exp_util./(1+exp_util);
    exp_util = exp(-(const+FV));         % this is now the utility of the external good
    prob=1./(1+exp_util);
    pmat = [prob 1-prob];
    pmat = pmat.*choice_dv;
    [r c p] = find(pmat);                                             % I*1
    LL = -sum(log(p));                                                % 1*1
    
    LL_sum = LL_sum + LL;
    
end

% format long g
% LL

%  p_all = utility_all./repmat(sum(utility_all,2),1,J);
%  d = p_all-choice_dv;                     % beware the sequence of what minus what !!
%  
%  Gt = zeros(I,(2*(J-1)+1));                         % number of variables
% %  Gt(:,1) = sum(ones(I,1).*d(:,1), 2);
% %  Gt(:,2) = sum(ones(I,1).*d(:,2), 2);
% %  Gt(:,3) = sum(ones(I,1).*d(:,3), 2);
% %  Gt(:,4) = sum(ones(I,1).*d(:,4), 2);
% for j = 1:(J-1)
%     Gt(:,j) = sum(ones(I,1).*d(:,j), 2);
%     Gt(:,j+1+J-1) = sum(demograph.*d(:,j), 2);
% end
% Gt(:,5) = sum(IV2.*d, 2);
% %  Gt(:,6) = sum(demograph.*d(:,1), 2);
% %  Gt(:,7) = sum(demograph.*d(:,2), 2);
% %  Gt(:,8) = sum(demograph.*d(:,3), 2);
% %  Gt(:,9) = sum(demograph.*d(:,4), 2);
%  
%  gr = sum(Gt,1);
%  
%  H = inv(Gt'*Gt);                          % H is actually the inverse of B here (B is approaching -Ht asymptotically) 
%  disp(gr);
%  se = sqrt(diag(H));                      % H is the covariance matrix here ? 
%  % disp(se);
%  % fprintf('standard error %d ',se);

end