% w: not working under MATLAB 2016      
% Generate Weibull data with A depending on predictor X
      x = 4*rand(100,1); A = 50*exp(-0.5*x); B = 2;
      y = wblrnd(A,B);
   
      % Fit Cox model
      [b,logL,H,st] = coxphfit(x,y);
   
      % Show Cox estimate of baseline survivor and known Weibull function
      figure;
      stairs(H(:,1),exp(-H(:,2)))
      xx = linspace(0,100);
      line(xx,1-wblcdf(xx,50*exp(-0.5*mean(x)),B),'color','r')
      title(sprintf('Baseline survivor function for X=%g',mean(x)));

      % Assess model adequacy by residual plots
      load(fullfile(matlabroot,'examples','stats','readmissiontimes.mat'))
      [b1,logL1,H1,st1] = coxphfit([Age Sex Weight Smoker],ReadmissionTime,...
                          'censoring',Censored,'ties','efron');

      % Draw a deviance residual plot
      Deviance = st1.devres;
      SmokerDev = Deviance(Smoker==1);
      NonsmokerDev = Deviance(Smoker==0);
      figure;
      plot((1:length(SmokerDev)),SmokerDev,'ro',...
           (1:length(NonsmokerDev)),NonsmokerDev,'b*')
      % Add a horizontal reference line with intercept 0
      hline = refline([0,0]);
      hline.Color = 'k';
      legend('Smoker','Non-smoker','Location','Best')
      xlabel('Observation number')
      ylabel('Deviance residuals')

      % Draw a martingale residual plot
      Schres = st1.schres(:,4);      
      figure;
      plot(ReadmissionTime(Smoker==1),Schres(Smoker==1),'ro',...
           ReadmissionTime(Smoker==0),Schres(Smoker==0),'b*')
      hline = refline([0,0]);
      hline.Color = 'k';
      legend('Smoker','Non-smoker','Location','Best')
      xlabel('Time')
      ylabel('Schoenfeld residuals')
      
      % Refit the model by stratifying a covariate
      [b2,logL2,H2,st2] = coxphfit([Age Sex Weight], ReadmissionTime,...
                          'censoring',Censored,'ties','efron','strata',Smoker);