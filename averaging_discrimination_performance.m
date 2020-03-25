%% averaging performance across blocks of trials
p=1;
for j = 23:28
    if j == 25 %|| j == 26
        continue
    end
    for i=1:5
        correct(:,i) = vertcat(total_perf(j).trial((20*i)-19:20*(i)).correct_choice);
        
    end
    avg(p,:) = sum(correct)/length(correct)
    p=p+1;
    clear correct
end
z(4,:) = mean(avg)

bar(mean(z))
ylim([.5 1])
title('Performance across 20 trial blocks')
xlabel('Block')
ylabel('Performance')
bins = 20*(1:5)
correct = vertcat(total_perf(28).trial(bins).correct_choice)
subarray(
excel_tastes = {'Suc','NaCl'};
excel_directions = [1 2];

 file = dir('*.rhd');
 [~,~, ~,~,summary]=process_intan_v4_behavior_only(file.name,excel_tastes, excel_directions);
 
 notfirst = vertcat(total_perf(27).trial(50:end-20).correct_choice);
 nfavg = sum(notfirst)/length(notfirst)
 %% moving average
 j = 24;
 p=1;
 numtrialsperblock=20;
 figure;
 for j=11-5:11
     for i =1:length(total_perf(j).trial)-numtrialsperblock
         moving_avg(i) = sum(vertcat(total_perf(j).trial((i:numtrialsperblock-1+i)).correct_choice))/numtrialsperblock;
     end
     hold on
     subplot(3,2,p);
     plot(moving_avg)
     p=p+1;
     clear moving_avg

title('Moving average')
xlabel('Blocks')
ylabel('Performance')
ylim([.5 1])
 end
 suptitle(append(total_perf(1).mouseID,' -- Pure Tastes'))
 
 %% avg performance
 j = 24;
 p=1;
 numtrialsperblock=20;
 figure;
 for j=length(total_perf)-10:length(total_perf)
   avg(p) = total_perf(j).total_performance
   bias(p)=total_perf(j).bias
   p=p+1;


 end
 
 subplot(2,1,1);
 scatter(1:length(avg),avg)
 title('Average Performance')
ylabel('Performance')
ylim([.3 1])
line([0 length(avg)],[.5 .5],'LineStyle','--')
line([0 length(avg)],[.7 .7],'Color','r','LineStyle','--')
 subplot(2,1,2);
 scatter(1:length(bias),bias)
 title('Bias')
 xlabel('Session')
ylabel('Bias')
ylim([-1 1])
line([0 length(bias)],[0 0],'LineStyle','--')
 suptitle('JK092 -- Mixtures')
%% stats
n=1:130;
p=zeros(1,130);
p(:)=.5;
simulated_trials = zeros(1000,130);
for i=1:1000
simulated_trials(i,:) = binornd(1,p);
end
for j = 1:1000
for i =1:length(simulated_trials(1,:))-20
   moving_avg(j,i) =  sum(simulated_trials(j,i:19+i))/20;
end
end
rand_avg = mean(moving_avg);
rand_stdev = 2*std(moving_avg);
rand_sem = 2*std(moving_avg)./1000;
plot(rand_avg)
errorbar(rand_avg,rand_stdev)
ylim([0 1])
title('Random blocks of 20 trials, error bars  are 2 standard deviations')
xlabel('Trials')
ylabel('Performance')
 %% glm testing
 n=28;
 prev_dir = [];
 prev_sucess=[];
 blocks=[];
 bins=[];
 prev_sucess = vertcat(total_perf(n).trial(1:end-1).correct_choice)';
 
 for i =2:length(total_perf(n).trial)
     if total_perf(n).trial(i).L_R_trial == total_perf(n).trial(i-1).L_R_trial
    prev_dir(1,i-1) = 1;
     else
    prev_dir(1,i-1) = 0;
     end
 end
 
 num_block = round(length(total_perf(n).trial)/20)+1;
 bins = (1:num_block)*20;
 blocks = zeros(num_block,length(total_perf(n).trial));
 for i =1:length(bins)
    
    blocks(i,(20*i)-19:20*i) =  1;
 end
 blocks(:,length(total_perf(n).trial)+1:end) = [];
 blocks(:,1) = [];
 
 
 x = vertcat(prev_dir, prev_sucess,blocks)';
 y = vertcat(total_perf(n).trial(2:end).correct_choice);
 
 [b,dev,stats] = glmfit(x,y,'binomial');
%  yfit = glmval(b,X,'probit','size');
% plot(x(:,2), y,'o',x(:,2),yfit,'-','LineWidth',2)