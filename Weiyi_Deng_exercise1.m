% Weiyi Deng 401716
c = 1:100;
i=1;
for x = 2:100
if isprime(x)==1
y = x*c;
count = sum(y<=100);
total_count(i)= count;         % the ith element of total_count is the current count
else
end
end
disp(total_count)
bar(total_count)
num_prime = sum(isprime(2:100));
set(gca,'XTick',1:num_prime)
set(gca, 'xticklabel', {primes(100)})
title('Primes within 100')
xlabel('Primes')
ylabel('Counts being a factor')


        