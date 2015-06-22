n = 2:100
%primes = []; 
for i = 1:length(n)
    if n(i)~=0          % n(i) stands for the ith element in vector n !!
    for j = i + 1:length(n)
        if mod(n(j),n(i))==0
            n(j) = 0
        end
    end
    end
    if n(i)~=0
        primes(length(primes)+1) = n(i)
    end
end
primes