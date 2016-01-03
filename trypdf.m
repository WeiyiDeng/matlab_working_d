X = -2.5:0.01:2.5;
Y = normpdf(X,0,1);
plot(X,Y)

x_ea = norminv(0.16,0,1);
x_ma = norminv(0.84,0,1);

line([x_ea x_ea],[0 0.5]);
line([x_ma x_ma],[0 0.5]);
ylim([0 1])
line([-2.5 2.5],[0 0]);
xlabel([]);
ylabel([]);
set(gca,'Visible','off')