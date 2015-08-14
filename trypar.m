% clc
% clear all

% feature('numCores')

% w: turn on matlabpool to use parfor with multiple cores
% matlabpool local 2               % same as parpool in new version matlab

num = 300
y = cell(num,1);
parfor i = 1:num    
    stuff = round(inv(rand(100)).*1000);
    b = [];
    for j = 1:1000
        bt = find(stuff==j);
        b = [b bt']; 
    end
    y{i} = b;
end

% matlabpool CLOSE

% job = batch('trypar')
% delete(job)                       % after job finishes



