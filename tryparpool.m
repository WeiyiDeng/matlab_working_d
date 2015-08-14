% feature('numCores')
matlabpool local 6               % same as parpool in new version matlab

tic
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
toc

matlabpool CLOSE

% job = batch('trypar')
% delete(job)                       % after job finishes



