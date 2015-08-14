% build a struct s with one field and three values
s(1).f1 = rand(3, 6);
s(2).f1 = magic(12);
s(3).f1 = ones(5, 10);

% build an anonymous function to get number of elements
% variable funky is the function handle here
funky = @(x) numel(x.f1);
funky(s(1))

% % or alternatively
% funky2 = @(x) mean(mean(x));
% funky2(s(1).f1)

% no need to store the anonymous function with a variable 
% can just call the anonymous function within the arrayfun function 
% (can also do this in other functions)
% arrayfun will output the results of each element in the input array (in
% this case s) from the given function (in this case @(x))
arrayfun(@(x) numel(x.f1), s)