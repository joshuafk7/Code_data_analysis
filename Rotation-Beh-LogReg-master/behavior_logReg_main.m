% Main 
% evaluates 4 models and plots Log likelihood, AIC and beta weights.
% Model 1 - stimulus, model 2 - stim + choice history, model 3 - stim +
% reward hist, model 4 - reward history + stim history + reward hist * stim
% hist.

stim = []; rewardHist = []; choiceHist = []; 

trialStrength = [1,0.75,0.55,0,0.45,0.25,0];
 
for i = 2:size(trial,2)
     
    stim(i-1) = trialStrength(trial(i).tasteStim);
    
    if isempty(trial(i-1).rewardOnset)
        rewardHist(i-1) = 0;
    else
        rewardHist(i-1) = 1;
    end
    
    choiceHist(i-1) = -1*trial(i-1).choice;
    
end
 
xVar_1 = stim';
xVar_2 = [stim', choiceHist',rewardHist', choiceHist'.*rewardHist'];
xVar_3 = [stim',choiceHist'];
xVar_4 = [stim',rewardHist'];

yVar = [trial.choice]==-1;
yVar = yVar(2:end);

yLog_1 = [sum(yVar(xVar_1==0)==1)/sum((xVar_1==0)), sum(yVar(xVar_1==0.25)==1)/sum((xVar_1==0.25)), ...
              sum(yVar(xVar_1==0.45)==1)/sum((xVar_1==0.45)), sum(yVar(xVar_1==0.55)==1)/sum((xVar_1==0.55)),...
              sum(yVar(xVar_1==0.75)==1)/sum((xVar_1==0.75)), sum(yVar(xVar_1==1)==1)/sum((xVar_1==1))];


for mod = 1:4


eval(sprintf("[logitCoef_%d,dev_%d,stats_%d] = glmfit(xVar_%d,yVar','binomial','link','logit');",mod,mod,mod,mod));
eval(sprintf("yhat_%d = glmval(logitCoef_%d,xVar_%d,'logit');",mod,mod,mod));

eval(sprintf('yLog_avg_%d = [];',mod));

    for s = [7 6 5 3 2 1]
        
        eval(sprintf("yLog_avg_%d = [yLog_avg_%d, nanmean(yhat_%d(xVar_1==trialStrength(s)))];",mod,mod,mod));
                  
    end


end

figure(1)         
plot([0 0.25 0.45 0.55 0.75 1],yLog_1,'ko-','LineWidth',1.2);    

hold on;

plot([0 0.25 0.45 0.55 0.75 1],yLog_avg_1,'ro-',[0 0.25 0.45 0.55 0.75 1],yLog_avg_3,'bo-',...
    [0 0.25 0.45 0.55 0.75 1],yLog_avg_4,'yo-',[0 0.25 0.45 0.55 0.75 1],yLog_avg_2,'go-','LineWidth',1.2)
    
xlabel('Sucrose Concentration','FontSize',12)
ylabel('Choice Probability (Sucrose)','FontSize',12)

legend({'Data','Model 1','Model 2','Model 3','Model 4'},'location','SouthEast')
    

% AIC is computed as (-Log Likelihood + 2*(parameters))
aic = [dev_1+2, dev_3 + 4,dev_4 + 4, dev_2+8];

figure(2)
subplot(1,2,1)
bar([-dev_1, -dev_3, -dev_4,-dev_2])
ylabel('Log Likelihood Value','FontSize',12)
xticklabels({'Model 1','Model 2','Model 3','Model 4'})
ylim([-70 -50])

axis square;

subplot(1,2,2)
bar(aic)
ylabel('AIC','FontSize',12)
xticklabels({'Model 1','Model 2','Model 3','Model 4'})
ylim([55 70])

axis square;



figure(3)
bar(logitCoef_2(2:end))
ylabel('\beta Value','FontSize',12)
xticklabels({'Stimulus','Choice Hist','Reward Hist','Choice * Reward Hist'})
axis square;
