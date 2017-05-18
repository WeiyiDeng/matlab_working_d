% w: not working under MATLAB 2016a
%% try example 
% from https://www.mathworks.com/help/stats/cox-proportional-hazards-model-with-time-dependent-covariates.html
% prep data
A = [1       0        46      0            1       0.3
    2       0        50      1            0       0.2
    2      50       100      1            0      0.23
    2     100       138      1            0      0.39
    3       0        50      1            1      0.18
    3      50        94      0            1      0.22
    4       0        50      1            0      0.21
    4      50        50      0            0       0.2
    5       0        50      1            0      0.25
    5      50       100      1            0      0.21
    5     100       106      0            0      0.42
    6       0        50      1            0      0.21
    6      50        98      0            0      0.22];

ID = A(:,1);
tStart = A(:,2);
tStop = A(:,3);
Censoring = A(:,4);
Sex = A(:,5);
Lab = A(:,6);

labCP = table(ID,tStart,tStop,Censoring,Sex,Lab)

% 
idxInvalid = labCP.ID(find(labCP.tStart == labCP.tStop))

labCP(find(labCP.ID==idxInvalid),:)

idxAdjust = find(labCP.ID==idxInvalid);
labCP.tStop(idxAdjust(1)) = labCP.tStop(idxAdjust(1))-0.5;
labCP.tStart(idxAdjust(2)) = labCP.tStart(idxAdjust(2))-0.5;
labCP(idxAdjust,:)

% 
X = [labCP.Sex labCP.Lab];
T = [labCP.tStart labCP.tStop];
b = coxphfit(X,T,'Censoring',labCP.Censoring,'Baseline',0)
% gives error message: Error using coxphfit (line 70) Y must be a real vector.
% try in 2017a ?