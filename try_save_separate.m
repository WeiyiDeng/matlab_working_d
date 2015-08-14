for i = 1:10
    A = randn(5);
    save(sprintf('trysavefiles_%02d',i), 'A');
end
    
    
    
% for i = 1:10
% 
% a = i;
% 
% saving_variable_name = sprintf('Results %d',i); 
% save(saving_variable_name);
% 
% end