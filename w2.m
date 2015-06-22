x = linspace(-2.1, 0.6, 50);
y = linspace(-1.1, 1.1, 50);
for i = 1:length(x)
    for j = 1:length(y)
        c = complex(x(i),y(j));
        z = 0;
        n = 1;
        [~, num_itera] = itera(z,c,n);
        num_itera
        combine(j,i)=num_itera;
    end
end
imagesc(combine)
