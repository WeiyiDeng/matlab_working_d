%// example data
clear all
rng(0,'twister')
data = randn(1000,1);

% w: prevent the axes from being skewed
x = linspace(-4,4,100);
% x_q = quantile(x,0.2:0.1:0.8);
% y = 16 - x_q.^2;
% y = 16 - x.^2;
y = -2.*ones(size(x));

%// generate two axes at same position
ax1 = axes;
ax2 = axes('Position', get(ax1, 'Position'),'Color','none');

%// move second axis to the right, remove x-ticks and labels
set(ax2,'YAxisLocation','right')
set(ax2,'XTick',[])

%// plot hist and line plot
hist(ax1,data,50); hold on
% plot(ax2,x_q,y)
plot(ax2,x,y,'k')

ylabel(ax1,'label of hist')
ylabel(ax2,'label of plot')
xlabel(ax1,'Hello World!')

x_point = quantile(data(:,1),[0.16 0.84]);                           %<- to change
% x_point = x_point +[0.06, -0.06];
% line([x_point(1) x_point(1)],[0 16*10^5], 'Color', 'r','Linestyle',':')
% line([x_point(2) x_point(2)],[0 16*10^5], 'Color', 'r','Linestyle','--')
line([x_point(1) x_point(1)],[-2 16], 'Color', 'r','Linestyle','--')
line([x_point(2) x_point(2)],[-2 16], 'Color', 'r','Linestyle','--')
text(x_point(1)+0.1,0.001,'16%')
text(x_point(2)+0.1,0.001,'84%')

x_q = quantile(x,0.2:0.1:0.5);                                % w: now can change the quantile as wish without the lines being shifted
y = 16 - x_q.^2;
plot(ax2,x_q,y)