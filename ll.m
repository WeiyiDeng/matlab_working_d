function LL = ll(b)

heat = csvread('data5.csv');
I = 900;
J = 5;

choicemat = repmat(heat(:,2),1,J);
testmat = repmat(1:J,I,1);
choice_dv = choicemat==testmat;

ic = heat(:,3:7)./100;
oc = heat(:,8:12)./100;

utility_all = exp(ic.*b(1)+oc.*b(2));
pmat = utility_all.*choice_dv./repmat(sum(utility_all,2),1,J);
[r c p] = find(pmat);
LL = -sum(log(p));

format long g
LL

end