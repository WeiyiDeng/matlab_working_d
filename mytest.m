t = 0:0.01:3

for k = 0:0.1:1
    y = 1./exp(k).^t;
    plot(t,y)
    hold on
end
hold off

plot(t,1./5.^t)
hold on
plot(t,1./0.3.^t)
hold off