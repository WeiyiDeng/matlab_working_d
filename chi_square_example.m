% chi-square example with parameter k = 1,2,3,4
x = 0:0.2:15;
y = chi2pdf(x,4);
plot(x,y)
hold on
y = chi2pdf(x,3);
plot(x,y)
y = chi2pdf(x,2);
plot(x,y)
y = chi2pdf(x,1);
plot(x,y)
hold off

x = 0:0.2:15;
y = gampdf(x,2,0);
plot(x,y)
hold off

x = 0:0.2:30;
for k = 0:10
    for theta = 1:5
        y = gampdf(x,k,theta);
        plot(x,y)
        hold on
    end
end
hold off
    