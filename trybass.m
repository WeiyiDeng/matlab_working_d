% w: can also use excel spreadsheet

m = 30;
p = 0.0112;
q = 0.5;
% p = 0.5;
% q = 0.001;
% p = 0.1;
% q = 0.1;

T = 24;
N = zeros(T,1);
S = zeros(T,1);

N(1) = m*p;
S(1) = N(1);

for t=2:T
    S(t) = (p + q/m*N(t-1))*(m-N(t-1));
    N(t) = N(t-1)+S(t);
end

plot(S)
% hold on
% plot(N)
% hold off