% clear all
% 
% tic
% test = cell(4,1);
% for N = 1:4
%     % createTask(job, @task, #outputs, {inputs(N)})
%     test{N} = tryparfun(N*100, 1000);
% end
% toc


job = batch('trypar','matlabpool',1)

wait(job)
load(job)
% load(job,'A')
% plot(A)

delete(job)
% destroy(job)
clear job

