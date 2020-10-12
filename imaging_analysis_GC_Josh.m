% Imaging analysis for GC project
%% rigister your image
 %edit regtiff
%% load imaging data and ROI

% neuron = results;
 
load('neuron_filtered.mat','neuron')
neuron = neuron_new;
clearvars -except neuron
%% load the event
% file = dir('*.rhd');
% a = {'Suc','S_75_25','S_55_45','S_45_55','S_25_75', 'NaCl'};
% % a = {'Suc','S_85_15','S_65_35','S_35_65','S_15_85', 'NaCl'};
% b = [2 2 2 1 1 1];
% [tastes,~,~,trial,summary] = process_intan_v4_behavior_only(file, a, b);
% clear file a b;
fldr = pwd;
date=fldr(end-5:end);
date_2=str2num(date);
cd ..
cd summary
load('total_perf.mat','total_perf')
   dates = zeros(length(total_perf));
   dates = vertcat(total_perf.date);
   dates = str2num(dates);
   
current_session=find(dates==date_2);
tastes=total_perf(current_session).tastes;
trial=total_perf(current_session).trial;
summary=total_perf(current_session);
cd ..

cd(date)
clear current_session date date_2  dates fldr
%% extract traces and spikes to trial structurebased on index of frames in C from CNMF   
for i =1:length(trial)
    trial(i).traces = neuron.C(:,trial(i).imaging_frames_index);
    trial(i).spikes = full(neuron.S(:,trial(i).imaging_frames_index));
end

%% reshape trials to neurons

neurons = trial2neuron5tastant_Josh(trial,tastes);

%%
% save('imaging_update','neurons','neuron','trial','tastes','summary')
%%
%% binned data - create 4d matrix 
%4D array called binnedC, dim1=trial, dim2=bin, dim3=taste, dim4=neuron
bins = -4:2:6; %bin size for comparison 2 second
binnedC = fun_binning(1,trial, bins, tastes, neurons);
plotbins = -4:.5:6; %smaller bin sizes for plotting
plotbinnedC = fun_binning(1,trial, plotbins, tastes, neurons);
%% binnedC2 - trial x bin x neuron
binnedC2 = fun_binning(2, trial, bins, tastes, neurons);
plotbinnedC2 = fun_binning(2,trial, plotbins, tastes, neurons);
%% binnedC3 trial x bin x correct vs error (L/R) x neuron
binnedC3 = fun_binning(3, trial, bins, tastes, neurons);
plotbinnedC3 = fun_binning(3,trial, plotbins, tastes, neurons);
%% binnedC4 trial x bin x taste, correct vs error x neuron
binnedC4 = fun_binning(4, trial, bins, tastes, neurons);
plotbinnedC4 = fun_binning(4,trial, plotbins, tastes, neurons);
%% binnedC5 trial x bin x L/R x neuron
binnedC5 = fun_binning(5, trial, bins, tastes, neurons);
plotbinnedC5 = fun_binning(5,trial, plotbins, tastes, neurons);

%% binnedC6 trial x bin x taste(correct_only) x neuron
binnedC6 = binnedC4(:,:,1:length(tastes),:);
plotbinnedC6 = plotbinnedC4(:,:,1:length(tastes),:);

%% binnedC7 trial x bin x correct(2), error (1) x neuron
binnedC7 = fun_binning(7, trial, bins, tastes, neurons);
plotbinnedC7 = fun_binning(7,trial, plotbins, tastes, neurons);

%% check multiple time bins against baseline for taste responsiveness
baseline =2;
comparison = 3;
taste_resp_neurons = fun_taste_responses(neurons, tastes, binnedC6, binnedC7(:,:,2,:), baseline, comparison);
% plot_sig_neurons(taste_resp_neurons, plotbins, plotbinnedC2,summary)

comparison = 4;
delay_resp_neurons = fun_taste_responses(neurons, tastes, binnedC6, binnedC7(:,:,2,:), baseline, comparison);
% plot_sig_neurons(taste_resp_neurons, plotbins, plotbinnedC2,summary)

comparison = 5;
choice_resp_neurons = fun_taste_responses(neurons, tastes, binnedC6, binnedC7(:,:,2,:), baseline, comparison);
% plot_sig_neurons(choice_resp_neurons, plotbins, plotbinnedC2,summary)
%% code for SVM classification
test=[];
test2=[];
x=[100 75 55 45 25 0;1 2 3 4 5 6];
% test=zscore(test);
for i=1:length(trial)
    for j=1:length(tastes)
     if ~isempty(find(convertCharsToStrings(trial(i).TasteID)==tastes(j)))
         trial_labels(1,i)=x(2,j);
     end
    end
end

[numtastes(:,2),numtastes(:,1)]=groupcounts(test2');
subsample_num_trials=min(numtastes(:,2));

randomization=randperm(length(trial));
for i=1:length(trial)
    shuff_trials(i)=trial(randomization(i));
end

%% subsample so there are equal number of trials per taste
trials_deletion=[];
taste_counts=zeros(1,length(tastes));
p=1;
for i=1:length(shuff_trials)
    for j=1:length(tastes)
        if ~isempty(find(convertCharsToStrings(shuff_trials(i).TasteID)==tastes(j)))
            %          test2(1,i)=x(2,j);
            taste_counts(j)=taste_counts(j)+1;
              if taste_counts(j)>subsample_num_trials
                trials_deletion(p)=i;
                p=p+1;
            end
            
        end
      
    end
end
randomization(trials_deletion)=[];

for i=1:length(shuff_trials)
    binnedC2_subsample(i,:) = binnedC2(randomization(i),3,:);
end

shuff_trials(trials_deletion)=[];
%
test2=[];
for i=1:length(shuff_trials)
    for j=1:length(tastes)
     if ~isempty(find(convertCharsToStrings(shuff_trials(i).TasteID)==tastes(j)))
         test2(1,i)=x(1,j);
     end
    end
end
test2=test2';
test1=[];
test1=zscore(binnedC2_subsample(:,clust1));
% [numtastes(:,2),numtastes(:,1)]=groupcounts(test2')
% test2=test2';


%% population averages - zscorred, correct only different epochs
figure(789)

subplot(4,1,1)
popavgtaste1 = mean(zscore(squeeze(nanmean(plotbinnedC7(:,:,2,taste_resp_neurons),1)),0,1),2)';
hold on
plot(plotbins(1:end-1),popavgtaste1);
y=ylim;
rectangle('Position',[0 y(1) 2 y(2)-y(1)],'EdgeColor', 'none', 'FaceColor', [.8, .8, .8, 0.5])
line([0 0],[y(1) y(2)],'Color','r','LineStyle','--')
title('Taste Responsive Neurons')
% xlabel('Time (sec)')
ylabel('Spiking Z score')

subplot(4,1,2)
y=[];
popavgtaste2 = mean(zscore(squeeze(nanmean(plotbinnedC7(:,:,2,delay_resp_neurons),1)),0,1),2)';
hold on
plot(plotbins(1:end-1),popavgtaste2);
y=ylim;
rectangle('Position',[2 y(1) 2 y(2)-y(1)],'EdgeColor', 'none', 'FaceColor', [.8, .8, .8, 0.5])
line([0 0],[y(1) y(2)],'Color','r','LineStyle','--')
title('Delay Responsive Neurons')
% xlabel('Time (sec)')
ylabel('Spiking Z score')

subplot(4,1,3)
y=[];
popavgtaste3 = mean(zscore(squeeze(nanmean(plotbinnedC7(:,:,2,choice_resp_neurons),1)),0,1),2)';
hold on
plot(plotbins(1:end-1),popavgtaste3);
y=ylim;
rectangle('Position',[4 y(1) 2 y(2)-y(1)],'EdgeColor', 'none', 'FaceColor', [.8, .8, .8, 0.5])
line([0 0],[y(1) y(2)],'Color','r','LineStyle','--')
title('Choice Responsive Neurons')
% xlabel('Time (sec)')
ylabel('Spiking Z score')

subplot(4,1,4)
plot_licking(summary)
set(gcf, 'Renderer', 'painters');

%%
interneurons_threshold = unique(overlapping_pair_filtered2(:,2));
% excitatory = 
% 
test555=intersect(interneurons_threshold,taste_resp_neurons);
test666=intersect(interneurons_threshold,delay_resp_neurons);
test777=intersect(interneurons_threshold,choice_resp_neurons);
task_related_interneurons_threshold = unique([test555;test666;test777]);

task_related_neurons_all = unique([taste_resp_neurons;delay_resp_neurons;choice_resp_neurons]);
task_related_excitatory = intersect(task_related_neurons_all, excitatory);


%%
figure(789)

subplot(2,1,1)

popavgtaste1 = mean(squeeze(nanmean(plotbinnedC7(:,:,2,interneurons_threshold),1)),2)';
popavgtaste2 = mean(squeeze(nanmean(plotbinnedC7(:,:,2,excitatory),1)),2)';
plot(plotbins(1:end-1),popavgtaste1);
% ylim([0.005 0.019])
hold on
plot(plotbins(1:end-1),popavgtaste2);
% ylim([0.005 0.019])
% ylim([0.009  0.016])
% ylim([0 0.05])
y=ylim;
rectangle('Position',[0 0 2 y(2)],'EdgeColor', 'none', 'FaceColor', [.8, .8, .8, 0.5])
line([0 0],[0 y(2)],'Color','r','LineStyle','--')
title('Population Activity of All Exh/Inh Neurons')
% xlabel('Time (sec)')
ylabel('Mean Spiking')
legend('Interneurons','All Neurons')



subplot(2,1,2)
plot_licking(summary)
set(gcf, 'Renderer', 'painters');
%%
figure(795)

subplot(2,1,1)

popavgtaste1 = mean(squeeze(nanmean(plotbinnedC7(:,:,2,task_related_interneurons_threshold(1:3)),1)),2)';
popavgtaste2 = mean(squeeze(nanmean(plotbinnedC7(:,:,2,task_related_excitatory),1)),2)';
plot(plotbins(1:end-1),popavgtaste1);
% ylim([0.005 0.019])
hold on
plot(plotbins(1:end-1),popavgtaste2);
% ylim([0.005 0.019])
% ylim([0.009  0.016])
% ylim([0 0.05])
y=ylim;
rectangle('Position',[0 0 2 y(2)],'EdgeColor', 'none', 'FaceColor', [.8, .8, .8, 0.5])
line([0 0],[0 y(2)],'Color','r','LineStyle','--')
title('Population Activity of Task Responsive Neurons')
% xlabel('Time (sec)')
ylabel('Mean Spiking')
legend('Interneurons','All Neurons')



subplot(2,1,2)
plot_licking(summary)
set(gcf, 'Renderer', 'painters');

%% plotting interneuron activity
% figure(7875)

% subplot(2,1,1)
popavgtaste1 = mean(squeeze(nanmean(binnedC7(:,:,2,task_related_interneurons_threshold(1:3)),1)),2)'; %correct
popavgtaste2 = mean(squeeze(nanmean(binnedC7(:,:,1,task_related_interneurons_threshold(1:3)),1)),2)'; %error
popavgtaste3 = mean(squeeze(nanmean(binnedC7(:,:,2,task_related_excitatory),1)),2)'; %correct
popavgtaste4 = mean(squeeze(nanmean(binnedC7(:,:,1,task_related_excitatory),1)),2)'; %error

figure(884)
x={'baseline1','baseline2','sampling','delay','choice'};
y=[popavgtaste1;popavgtaste2;popavgtaste3;popavgtaste4]';
bar(categorical(x,x),y)
legend('Correct interneurons','Error interneurons', 'Correct Excitatory','Error Excitatory','Location','northwest')
ylabel('Inferred Spiking (a.u.)')
title('Mean Spiking - Correct vs Error Trials')
hold on
%% stats
comparison_bin=3;

tastebin_corr_interneurons = nanmean(squeeze(binnedC7(:,comparison_bin,2,task_related_interneurons_threshold(1:3))),2); %correct
tastebin_err_interneurons = nanmean(squeeze(binnedC7(:,comparison_bin,1,task_related_interneurons_threshold(1:3))),2); %error

ranksum(tastebin_corr_interneurons,tastebin_err_interneurons) %comparison of mean activity of task related interneuron between correct and error


tastebin_cor_allneurons = nanmean(squeeze(binnedC7(:,comparison_bin,2,task_related_neurons_all)),2); %correct
tastbin_err_allneurons = nanmean(squeeze(binnedC7(:,comparison_bin,1,task_related_neurons_all)),2); %error

ranksum(tastebin_cor_allneurons,tastbin_err_allneurons) %comparison of mean activity of all task related neurons between correct and error

%% looking at interneurons by taste correct vs error
comparison = 3;
avg_by_taste_corr = mean(squeeze(nanmean(binnedC4(:,comparison,1:6,task_related_interneurons_threshold(1:3)),1)),2)'; %error
avg_by_taste_err = mean(squeeze(nanmean(binnedC4(:,comparison,7:12,task_related_interneurons_threshold(1:3)),1)),2)'; %correct
x=[100 75 55 45 25 0];
y=[avg_by_taste_corr;avg_by_taste_err];
% close(643)
figure(643)
bar(categorical(x,x),y);
legend('correct','error')
ylabel('Inferred Spiking (a.u.)')
xlabel('Sucrose Concentration')
title('Interneuron Responses in Sampling Bin')
%% %% looking at interneurons by taste correct only
avg_by_taste_corr = mean(squeeze(nanmean(binnedC4(:,3,1:6,task_related_interneurons_threshold(1:3)),1)),2)'; %error
avg_by_taste_err = mean(squeeze(nanmean(binnedC4(:,3,7:12,task_related_interneurons_threshold(1:3)),1)),2)'; %correct
x=[100 75 55 45 25 0];
y=[avg_by_taste_corr;avg_by_taste_err];
close(643)
figure(643)
bar(categorical(x,x),y);
legend('correct','error')
ylabel('Inferred Spiking (a.u.)')
xlabel('Sucrose Concentration')
title('Interneuron Responses in Sampling Bin')
%% stats
comparison_bin=3;

Suc55_corr_interneurons = nanmean(squeeze(binnedC4(:,comparison_bin,3,task_related_interneurons_threshold(1:3))),2); %correct
Suc55_err_interneurons = nanmean(squeeze(binnedC4(:,comparison_bin,9,task_related_interneurons_threshold(1:3))),2); %error

ranksum(Suc55_corr_interneurons,Suc55_err_interneurons) %comparison of mean activity of task related interneuron between correct and error


% comparison_bin=3;

Suc55_corr_interneurons = nanmean(squeeze(binnedC4(:,comparison_bin,4,task_related_interneurons_threshold(1:3))),2); %correct
Suc55_err_interneurons = nanmean(squeeze(binnedC4(:,comparison_bin,10,task_related_interneurons_threshold(1:3))),2); %error

ranksum(Suc55_corr_interneurons,Suc55_err_interneurons) %comparison of mean activity of task related interneuron between correct and error
%% random subsampling of neurons vs interneurons
comparison_bin=3;
num_interneurons = length(task_related_interneurons_threshold);
p_shuff=[];
for i=1:1000
    subset=[];
    subset = randperm(length(task_related_excitatory),num_interneurons);
    tastebin_corr_random = nanmean(squeeze(binnedC7(:,comparison_bin,2,task_related_excitatory(subset))),2); %correct
    tastebin_err_random = nanmean(squeeze(binnedC7(:,comparison_bin,1,task_related_excitatory(subset))),2); %error

    p_shuff(i) = ranksum(tastebin_corr_random,tastebin_err_random); %comparison of mean activity of task related interneuron between correct and error

    tastebin_corr_rand_avg = nanmean(nanmean(squeeze(binnedC7(:,comparison_bin,2,task_related_excitatory(subset))),2)); %correct
    tastebin_err_rand_avg = nanmean(nanmean(squeeze(binnedC7(:,comparison_bin,1,task_related_excitatory(subset))),2)); %error
    difference(i) = tastebin_corr_rand_avg-tastebin_err_rand_avg;
end
find(p_shuff<0.02)

%%
comparison_bin=3;

    tastebin_corr_interneurons_avg = nanmean(nanmean(squeeze(binnedC7(:,comparison_bin,2,task_related_interneurons_threshold)),2)); %correct
    tastebin_err_interneurons_avg = nanmean(nanmean(squeeze(binnedC7(:,comparison_bin,1,task_related_interneurons_threshold)),2)); %error
    observed_difference=  tastebin_corr_interneurons_avg - tastebin_err_interneurons_avg;


stdev_x2 = std(difference)*2;
figure(6251)
histogram(difference,'Normalization','probability')
sdev=line([stdev_x2 stdev_x2],[0 .14],'Color','r','LineStyle','--');
% line([-stdev_x2 -stdev_x2],[0 .14],'Color','r','LineStyle','--')
obs = line([observed_difference observed_difference],[0 .14],'Color','g','LineStyle','--');
ylabel('Probability')
xlabel('Spiking Difference - Correct vs Error')
legend([sdev obs],'2 standard deviations','observed difference','Location','northwest')
title('1000 subsamples of task responsive neurons')
set(gcf,'Renderer','painters')
%%
% figure(7875)

% subplot(2,1,1)
popavgtaste1 = mean(squeeze(nanmean(binnedC7(:,:,2,test666),1)),2)'; %correct
popavgtaste2 = mean(squeeze(nanmean(binnedC7(:,:,1,test666),1)),2)'; %error
popavgtaste3 = mean(squeeze(nanmean(binnedC7(:,:,2,delay_resp_neurons),1)),2)'; %correct
popavgtaste4 = mean(squeeze(nanmean(binnedC7(:,:,1,delay_resp_neurons),1)),2)'; %error

figure(887)
x={'baseline1','baseline2','sampling','delay','choice'};
y=[popavgtaste1;popavgtaste2;popavgtaste3;popavgtaste4]';
bar(categorical(x,x),y)
legend('Correct interneurons','Error interneurons', 'Correct all neurons','Error all neurons')
ylabel('Inferred Spiking (a.u.)')
title('Mean Spiking Activity')
hold on
%%
% figure(7875)

% subplot(2,1,1)
popavgtaste1 = mean(squeeze(nanmean(binnedC7(:,:,2,test555),1)),2)'; %correct
popavgtaste2 = mean(squeeze(nanmean(binnedC7(:,:,1,test555),1)),2)'; %error
popavgtaste3 = mean(squeeze(nanmean(binnedC7(:,:,2,taste_resp_neurons),1)),2)'; %correct
popavgtaste4 = mean(squeeze(nanmean(binnedC7(:,:,1,taste_resp_neurons),1)),2)'; %error

figure(889)
x={'baseline1','baseline2','sampling','delay','choice'};
y=[popavgtaste1;popavgtaste2;popavgtaste3;popavgtaste4]';
bar(categorical(x,x),y)
legend('Correct interneurons','Error interneurons', 'Correct all neurons','Error all neurons')
ylabel('Inferred Spiking (a.u.)')
title('Mean Spiking Activity')
hold on
%%
% 
% figure(7873)
% 
% subplot(2,1,1)
% popavgtaste1 = mean(squeeze(nanmean(plotbinnedC7(:,:,2,test555),1)),2)';
% popavgtaste2 = mean(squeeze(nanmean(plotbinnedC7(:,:,1,test555),1)),2)';
% hold on
% plot(plotbins(1:end-1),popavgtaste1);
% plot(plotbins(1:end-1),popavgtaste2);
% hold off
% y=ylim;
% rectangle('Position',[0 y(1) 2 y(2)-y(1)],'EdgeColor', 'none', 'FaceColor', [.8, .8, .8, 0.5])
% line([0 0],[y(1) y(2)],'Color','r','LineStyle','--')
% title('Taste Responsive Neurons')
% xlabel('Time (sec)')
% subplot(2,1,2)
% plot_licking(summary)
% set(gcf, 'Renderer', 'painters');
% ylabel('Spiking Z score')


figure(7874)

subplot(2,1,1)
popavgtaste3 = mean(squeeze(nanmean(plotbinnedC7(:,:,2,test555),1)),2)'; %correct interneurons
popavgtaste5 = mean(squeeze(nanmean(plotbinnedC7(:,:,1,test555),1)),2)'; %correct all neurons

    
% popavgtaste4 = mean(squeeze(nanmean(plotbinnedC7(:,:,1,test666),1)),2)'; %error
hold on
plot(plotbins(1:end-1),popavgtaste3);
plot(plotbins(1:end-1),popavgtaste5);
% plot(plotbins(1:end-1),popavgtaste4);

y=ylim;
rectangle('Position',[2 0 2 y(2)],'EdgeColor', 'none', 'FaceColor', [.8, .8, .8, 0.5])
line([0 0],[y(1) y(2)],'Color','r','LineStyle','--')
title('All Interneurons')
xlabel('Time (sec)')
subplot(2,1,2)
plot_licking(summary)
set(gcf, 'Renderer', 'painters');


%%
plot(plotbins(1:end-1),popavgtaste1);
plot(plotbins(1:end-1),popavgtaste2);
plot(plotbins(1:end-1),popavgtaste3);
plot(plotbins(1:end-1),popavgtaste4);

y=ylim;
% rectangle('Position',[2 0 2 y(2)],'EdgeColor', 'none', 'FaceColor', [.8, .8, .8, 0.5])
line([0 0],[0 y(2)],'Color','r','LineStyle','--')
title('Delay Responsive Neurons')
xlabel('Time (sec)')
% subplot(2,1,2)
% plot_licking(summary)
set(gcf, 'Renderer', 'painters');
%% intersection analysis for venn diagram
 A= taste_resp_neurons;
 B=delay_resp_neurons;
 C=choice_resp_neurons;
 AB = intersect(A,B);
 AC = intersect(A,C);
 BC = intersect(B,C);
 ABC = intersect(AB,C);
 AandB=[A' B']';
 AandC=[A' C']';
 BandC=[B' C']';
 Aonly=setdiff(A, BandC);
 Bonly=setdiff(B, AandC);
 Conly=setdiff(C, AandB);
 venn(1) = length(Aonly);
 venn(2) = length(AB);
 venn(3) = length(Bonly);
 venn(4) = length(BC);
 venn(5) = length(Conly);
 venn(6) = length(AC);
 venn(7) = length(ABC);
 Allrespneurons=[];
 Allrespneurons = unique([A; B; C],'stable');
%  vennX(venn,.05)
 legend('A')
%%
 

 
%%
figure(221)

heatzscore1=zscore(squeeze(nanmean(plotbinnedC7(:,:,2,Allrespneurons),1))',0,2);
h1=heatmap(plotbins(1:end-1),Allrespneurons,heatzscore1);
h1.Colormap = parula;
h1.ColorLimits=[0 4];
title('Heatmap of Responsive Neurons');
xlabel('Time (sec)');
ylabel('Neurons');
set(gcf, 'Renderer', 'painters');

%%
 %% Clustering taste selective neurons - correct trials only
  p=[];
%   t = categorical({'Suc','S\_75\_25','S\_55\_45','S\_45\_55','S\_25\_75','NaCl'});
%  x=[100 75 55 45 25 0;];

  t = categorical({'S\_55\_45','S\_45\_55'});
 x=[ 55 45;];
  for i=1:length(taste_resp_neurons)
      taste_comparison=[];
  taste_comparison = squeeze(binnedC6(:,3,:,taste_resp_neurons(i)));
  [p,~,~] = anova1(taste_comparison,x,'off'); %compare all tastes in bin 3
  PP(i)=p;
  end
  idx1=find(PP<0.05);
  taste_selective=taste_resp_neurons(idx1);
  
  
%   intersect(taste_selective,interneurons)
%    taste_comparison=[];
%   taste_comparison = squeeze(binnedC(:,3,:,44));
%   [p,tbl,stats] = anova1(taste_comparison);
%   multcompare(stats)
%% 2way anova using smaller bins

% nan_values = find(isnan(plotbins_ANOVA));
% plotbins_ANOVA(nan_values)=[];
numbins=length(9:17);
tastelist=[];
p=1;
x=[100 75 55 45 25 0;1 2 3 4 5 6];
% test=zscore(test);
for i=1:length(trial)
    for j=1:length(tastes)
     if ~isempty(find(convertCharsToStrings(trial(i).TasteID)==tastes(j)))
         tastelist(1,p)=x(2,j);
         p=p+1;
     end
     
    end
end
tastelist=sort(tastelist);
tastelist=repmat(tastelist,1,numbins);

binlist = ones(1,length(tastelist));
p=1;
for j=1:numbins
    for i=1:length(trial)
        binlist(p) = j; 
        p=p+1;
    end
end


for q=1:length(taste_resp_neurons)
    plotbins_ANOVA=[];
    plotbins_ANOVA = reshape(squeeze(plotbinnedC6(:,9:17,:,taste_resp_neurons(q))),[],6);
    p=1;
for i=1:size(plotbins_ANOVA,1)
    for j=1:size(plotbins_ANOVA,2)
        if tastelist(i)==j
        plotbins_ANOVA2(p) = plotbins_ANOVA(i,j);
        p=p+1;
        end
    end
end



    
    
    
%     
% [numtastes(:,2),numtastes(:,1)]=groupcounts(test2');
% subsample_num_trials=min(numtastes(:,2));
% 
[pval44(:,q),~,~]=anovan(plotbins_ANOVA2,{tastelist, binlist},'display','off');
end
taste_resp_2wayANOVA=taste_resp_neurons(find(pval44(1,:)<0.05));

% plot_sig_neurons(taste_resp_2wayANOVA, plotbins, plotbinnedC7(:,:,2,:),summary)
%%

Q = squeeze(nanmean(binnedC6(:,3,:,taste_resp_neurons(idx1)),1));
Q=zscore(Q);
QQ = linkage(Q','Ward');
% % QQQ = cluster(QQ, 'Cutoff', 3, 'Depth',5);
QQQ = cluster(QQ, 'MaxClust', 4);

% neuron_number = 1:length(neurons);
clust1_neurons = taste_resp_neurons(idx1(find(QQQ == 1)));
clust2_neurons = taste_resp_neurons(idx1(find(QQQ == 2)));
clust3_neurons = taste_resp_neurons(idx1(find(QQQ == 3)));
clust4_neurons = taste_resp_neurons(idx1(find(QQQ == 4)));



% figure(899)
% dendrogram(QQ,'ColorThreshold',3,'Labels',num2str(taste_resp_neurons(idx1)),'Reorder',leaforder)
% title('Hierarchical Clustering of Taste Responsive Neurons')

% clust5_neurons = choice_resp_neurons(find(QQQ == 5));
% clust6_neurons = choice_resp_neurons(find(QQQ == 6));
% clust4_neurons = neuron_number(find(QQQ == 4));
% clust5_neurons = neuron_number(find(QQQ == 5));

%%

% colors=["b","r","g","c","m", "y","k"];
% figure(990)
% subplot(10,1,[1 2 3 4])
% dendrogram(QQ,'ColorThreshold',3,'Labels',num2str(taste_resp_neurons(idx1)),'Reorder',leaforder)
% set(gca,'fontsize',6)
% title('Hierarchical Clustering of Taste Responsive Neurons')
% clustplot1 = zscore(squeeze(nanmean(binnedC6(:,3,:,clustorder,1))));
% j=5;
% for i=1:6
% subplot(10,1,j)
% b(i)=bar(categorical(clustorder,clustorder),clustplot1(i,:),colors(i));
% ylabel(t(i))
% set(get(gca,'YLabel'),'Rotation',0)
% j=j+1;
% set(gca,'xtick',[])
% set(gca,'ytick',[])
% end
%% Plot normalized firing by cluster underneath dendrogram
figure(715)
D=pdist(Q');
leaforder=optimalleaforder(QQ,D);
temp1=taste_resp_neurons(idx1);
for i=1:length(temp1)
    clustorder(i) = temp1(leaforder(i));
end

subplot(10,1,[1 2 3 4])
dendrogram(QQ,'ColorThreshold',3,'Labels',num2str(taste_resp_neurons(idx1)),'Reorder',leaforder)
set(gca,'fontsize',6)
title('Hierarchical Clustering - Taste Selective Neurons')
clustplot1 = zscore(squeeze(nanmean(binnedC6(:,3,:,clustorder,1))));
p=5;
j=2;
i=1;
c=categorical(clustorder, clustorder);
% for j=1:6
%     subplot(10,1,p)
%     hold on

p=5;
for j=1:6
subplot(10,1,p)
    ylabel(t(j))
    
    set(get(gca,'YLabel'),'Rotation',0)
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    if j==6
       xlabel('Neurons') 
    end
for i=1:length(clustorder)
    hold on
   if intersect(clustorder(i), clust1_neurons) ~=0
     bar1= bar(c(i),clustplot1(j,i),'FaceColor',[0.4940 0.1840 0.5560]);
%       hold on
    clustupdate(i)=1;
   end
   if intersect(clustorder(i), clust2_neurons) ~=0
        hold on 
       bar1=bar(c(i),clustplot1(j,i),'c');
       clustupdate(i)=2;
%     ylabel(t(i))
%     set(get(gca,'YLabel'),'Rotation',0);
%    
%     set(gca,'xtick',[])
%     set(gca,'ytick',[])
   end
   if intersect(clustorder(i), clust3_neurons) ~=0
        hold on 
       bar1=bar(c(i),clustplot1(j,i),'r');
       clustupdate(i)=3;
%     ylabel(t(i))
%     set(get(gca,'YLabel'),'Rotation',0)
%     
%     set(gca,'xtick',[])
%     set(gca,'ytick',[])
   end
   if intersect(clustorder(i), clust4_neurons) ~=0
        hold on
       bar1=bar(c(i),clustplot1(j,i),'g');
       clustupdate(i)=4;
%     ylabel(t(i))
%     set(get(gca,'YLabel'),'Rotation',0)
%    
%     set(gca,'xtick',[])
%     set(gca,'ytick',[]) 
   end
end
p=p+1;
end
set(gcf, 'Renderer', 'painters');
%%
clusters = unique(clustupdate,'stable');
clust1=clustorder(find(clustupdate==clusters(1)));
clust2=clustorder(find(clustupdate==clusters(2)));
clust3=clustorder(find(clustupdate==clusters(3)));
clust4=clustorder(find(clustupdate==clusters(4)));

clustmean1=mean(zscore(squeeze(nanmean(binnedC6(:,3,:,clust1,1)))),2);
clustmean2=mean(zscore(squeeze(nanmean(binnedC6(:,3,:,clust2,1)))),2);
clustmean3=mean(zscore(squeeze(nanmean(binnedC6(:,3,:,clust3,1)))),2);
clustmean4=mean(zscore(squeeze(nanmean(binnedC6(:,3,:,clust4,1)))),2);
figure(666)
subplot(2,2,1)
scatter(x,clustmean1','g','filled')
title('Cluster 1')
xlabel('Sucrose Concentration')
ylabel('Spiking Z score')
subplot(2,2,2)
scatter(x,clustmean2','r','filled')
title('Cluster 2')
xlabel('Sucrose Concentration')
ylabel('Spiking Z score')
subplot(2,2,3)
scatter(x,clustmean3',[],[0.4940 0.1840 0.5560],'filled')
title('Cluster 3')
xlabel('Sucrose Concentration')
ylabel('Spiking Z score')
subplot(2,2,4)
scatter(x,clustmean4','c','filled')
title('Cluster 4')
xlabel('Sucrose Concentration')
ylabel('Spiking Z score')
% t = categorical({'Suc','S\_75\_25','S\_55\_45','S\_45\_55','S\_25\_75','NaCl'});
set(gcf, 'Renderer', 'painters');
%%
% figure(991)
% clustplot2 = squeeze(nanmean(binnedC6(:,3,:,clust2_neurons,1)));
% for i=1:6
% subplot(6,1,i)
% bar(categorical(clust2_neurons),clustplot2(i,:))
% end
% 
% figure(992)
% clustplot3 = squeeze(nanmean(binnedC6(:,3,:,clust3_neurons,1)));
% for i=1:6
% subplot(6,1,i)
% bar(categorical(clust3_neurons),clustplot3(i,:))
% end
% 
% figure(993)
% clustplot4 = squeeze(nanmean(binnedC6(:,3,:,clust4_neurons,1)));
% for i=1:6
% subplot(6,1,i)
% bar(categorical(clust4_neurons),clustplot4(i,:))
% end

%% 
figure(123);
subplot(4,1,1)
ylimit = plot_individual_neuron(1,Aonly(4), plotbins, plotbinnedC2);
hold on
rectangle('Position',[0 0 2 ylimit],'EdgeColor', 'none', 'FaceColor', [.8, .8, .8, 0.5])
hold off
subplot(4,1,2)
ylimit = plot_individual_neuron(1,Bonly(4), plotbins, plotbinnedC2);
hold on
rectangle('Position',[2 0 2 ylimit],'EdgeColor', 'none', 'FaceColor', [.8, .8, .8, 0.5])
hold off
subplot(4,1,3)
ylimit =plot_individual_neuron(1,Conly(5), plotbins, plotbinnedC2);
hold on
rectangle('Position',[4 0 2 ylimit],'EdgeColor', 'none', 'FaceColor', [.8, .8, .8, 0.5])
hold off
subplot(4,1,4)
plot_licking(summary)
set(gcf, 'Renderer', 'painters');
%%
figure(124);
subplot(3,1,1)
ylimit = plot_individual_neuron(1,112, plotbins, plotbinnedC2);
hold on
rectangle('Position',[0 0 2 ylimit],'EdgeColor', 'none', 'FaceColor', [.8, .8, .8, 0.5])
hold off
subplot(3,1,2)
ylimit = plot_individual_neuron(1,449, plotbins, plotbinnedC2);
hold on
rectangle('Position',[0 0 2 ylimit],'EdgeColor', 'none', 'FaceColor', [.8, .8, .8, 0.5])
hold off
subplot(3,1,3)
plot_licking(summary)
set(gcf, 'Renderer', 'painters');
 %% plot performance for psychometric
 
 figure
 psym=summary.ind_performance(1:3);
 psym(4:6) = 1-summary.ind_performance(4:6);
 scatter(x,psym)
 ylim([0 1])
 title('Session Performance')
 ylabel('Sucrose Choice')
 xlabel('Percent Sucrose')
 set(gcf, 'Renderer', 'painters');
 sigfunc = @(A,x) (A(1)./(1+exp(-A(2)*(A(3)+x))));
A_init = [1 0 1];
sigfunc(A_init,x)
A_fit = nlinfit(x, psym, sigfunc, A_init);
x_eval = (100:-5:0);
A_eval = sigfunc(A_fit,x_eval);
hold on
plot(x_eval,A_eval)
legend('Performance', 'Fit','Location','Northwest')
%% comparison between all taste pairs tastes- sampling
pval=[];
temp=[];
CC=[];
GC=[];
q=1;
% GC=nchoosek(x,2);
comparison=3;
[pval,tastepairs] = fun_compare_tastes(4,comparison, binnedC6, neurons, tastes);
GC=tastepairs;
[CC(:,1),CC(:,2)]=find(pval<0.05);
for i=1:size(CC,1)
    if ~isempty(intersect(CC(i,1),taste_resp_neurons))
        temp(q,1) = CC(i,1);
        temp(q,2) = CC(i,2);
        q=q+1;
    end
end
CC = temp;
selective_neurons = unique(CC(:,1));
numpairs = [1:size(tastepairs,1)];
[GG(:,2),GG(:,1)] = groupcounts(CC(:,2)) ;

paircomp=nan(6,6);
for i=1:length(GC)
paircomp(GC(i,1),GC(i,2)) = GC(i,3);
end
figure(8312)
h6=heatmap(t,t,paircomp');
h6.MissingDataColor = [1 1 1];
h6.MissingDataLabel='';
h6.GridVisible='off';
h6.ColorbarVisible='off';
 set(gcf, 'Renderer', 'painters');
%%
q=1;
pval=[];
temp=[];
CC=[];
[pval,tastepairs] = fun_compare_tastes(1,comparison, binnedC6, neurons, tastes);
[CC(:,1),CC(:,2)]=find(pval<0.05);
for i=1:size(CC,1)
    if ~isempty(intersect(CC(i,1),taste_resp_neurons))
        temp(q,1) = CC(i,1);
        temp(q,2) = CC(i,2);
        q=q+1;
    end
end
CC = temp;
plot_compare_tastes(comparison,plotbins, plotbinnedC6, CC)

%% clustering for neurons from 2way anova comparison
Q = squeeze(nanmean(binnedC6(:,3,:,taste_resp_2wayANOVA),1));
Q=zscore(Q);
QQ = linkage(Q','Ward');
% % QQQ = cluster(QQ, 'Cutoff', 3, 'Depth',5);
QQQ = cluster(QQ, 'MaxClust', 4);

% neuron_number = 1:length(neurons);
clust1_neurons = taste_resp_2wayANOVA(find(QQQ == 1));
clust2_neurons = taste_resp_2wayANOVA(find(QQQ == 2));
clust3_neurons = taste_resp_2wayANOVA(find(QQQ == 3));
clust4_neurons = taste_resp_2wayANOVA(find(QQQ == 4));
figure(716)
D=pdist(Q');
leaforder=optimalleaforder(QQ,D);
temp1=taste_resp_2wayANOVA;
for i=1:length(temp1)
    clustorder(i) = temp1(leaforder(i));
end

subplot(10,1,[1 2 3 4])
dendrogram(QQ,0,'ColorThreshold',7,'Labels',num2str(taste_resp_2wayANOVA),'Reorder',leaforder)
set(gca,'fontsize',6)
title('Hierarchical Clustering - Taste Selective Neurons')
clustplot1 = zscore(squeeze(nanmean(binnedC6(:,3,:,clustorder,1))));
p=5;
j=2;
i=1;
c=categorical(clustorder, clustorder);
% for j=1:6
%     subplot(10,1,p)
%     hold on

p=5;
for j=1:6
subplot(10,1,p)
    ylabel(t(j))
    
    set(get(gca,'YLabel'),'Rotation',0)
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    if j==6
       xlabel('Neurons') 
    end
for i=1:length(clustorder)
    hold on
   if intersect(clustorder(i), clust1_neurons) ~=0
     bar1= bar(c(i),clustplot1(j,i),'FaceColor',[0.4940 0.1840 0.5560]);
%       hold on
    clustupdate(i)=1;
   end
   if intersect(clustorder(i), clust2_neurons) ~=0
        hold on 
       bar1=bar(c(i),clustplot1(j,i),'g');
       clustupdate(i)=2;
%     ylabel(t(i))
%     set(get(gca,'YLabel'),'Rotation',0);
%    
%     set(gca,'xtick',[])
%     set(gca,'ytick',[])
   end
   if intersect(clustorder(i), clust3_neurons) ~=0
        hold on 
       bar1=bar(c(i),clustplot1(j,i),'c');
       clustupdate(i)=3;
%     ylabel(t(i))
%     set(get(gca,'YLabel'),'Rotation',0)
%     
%     set(gca,'xtick',[])
%     set(gca,'ytick',[])
   end
   if intersect(clustorder(i), clust4_neurons) ~=0
        hold on
       bar1=bar(c(i),clustplot1(j,i),'r');
       clustupdate(i)=4;
%     ylabel(t(i))
%     set(get(gca,'YLabel'),'Rotation',0)
%    
%     set(gca,'xtick',[])
%     set(gca,'ytick',[]) 
   end
end
p=p+1;
end
set(gcf, 'Renderer', 'painters');
%%
clusters = unique(clustupdate,'stable');
clust1=clustorder(find(clustupdate==clusters(1)));
clust2=clustorder(find(clustupdate==clusters(2)));
clust3=clustorder(find(clustupdate==clusters(3)));
clust4=clustorder(find(clustupdate==clusters(4)));

clustmean1=mean(zscore(squeeze(nanmean(binnedC6(:,3,:,clust1,1)))),2);
clustmean2=mean(zscore(squeeze(nanmean(binnedC6(:,3,:,clust2,1)))),2);
clustmean3=mean(zscore(squeeze(nanmean(binnedC6(:,3,:,clust3,1)))),2);
clustmean4=mean(zscore(squeeze(nanmean(binnedC6(:,3,:,clust4,1)))),2);
figure(667)
subplot(2,2,1)
scatter(x(1,:),clustmean1','g','filled')
title('Cluster 1')
xlabel('Sucrose Concentration')
ylabel('Spiking Z score')
subplot(2,2,2)
scatter(x(1,:),clustmean2','r','filled')
title('Cluster 2')
xlabel('Sucrose Concentration')
ylabel('Spiking Z score')
subplot(2,2,3)
scatter(x(1,:),clustmean3',[],[0.4940 0.1840 0.5560],'filled')
title('Cluster 3')
xlabel('Sucrose Concentration')
ylabel('Spiking Z score')
subplot(2,2,4)
scatter(x(1,:),clustmean4','c','filled')
title('Cluster 4')
xlabel('Sucrose Concentration')
ylabel('Spiking Z score')
% t = categorical({'Suc','S\_75\_25','S\_55\_45','S\_45\_55','S\_25\_75','NaCl'});
set(gcf, 'Renderer', 'painters');
%% comparison between tastes- sampling
pval=[];
temp=[];
CC=[];
comparison=3;
pval = fun_compare_tastes(1,comparison, binnedC6, taste_resp_neurons, tastes);
[CC(:,1),CC(:,2)]=find(pval<0.05);
q=1;
for i=1:size(CC,1)
    if ~isempty(intersect(CC(i,1),taste_resp_neurons))
        temp(q,1) = CC(i,1);
        temp(q,2) = CC(i,2);
        q=q+1;
    end
end
CC = temp;
plot_compare_tastes(comparison,plotbins, plotbinnedC6, CC)

%% compare left versus right trials
pval=[];
temp=[];
LR_Responsive=[];
comparison=3;
pval = fun_compare_tastes(3,comparison, binnedC, neurons, tastes);
LR_Responsive=find(pval<0.05);
q=1;
for i=1:size(LR_Responsive,1)
    if ~isempty(intersect(LR_Responsive(i,1),taste_resp_neurons))
        temp(q,1) = LR_Responsive(i,1);
%         temp(q,2) = LR_Responsive(i,2);
        q=q+1;
    end
end
LR_Responsive = temp;
%%
neuron2plot=733;
    figure
    hold on
    plot(plotbins(1:end-1), nanmean(plotbinnedC5(:,:,1,neuron2plot)))
    plot(plotbins(1:end-1), nanmean(plotbinnedC5(:,:,2,neuron2plot)))
    y=ylim;
    line([0 0],[0 y(2)],'Color','r','LineStyle','--')
    rectangle('Position',[0 0 2 y(2)],'EdgeColor', 'none', 'FaceColor', [.8, .8, .8, 0.5])
    title(append('Neuron ',num2str(neuron2plot,1)))
% plot_compare_tastes(comparison,plotbins, plotbinnedC, LR_Responsive)
%% compare between tastes in delay period
pval=[];
temp=[];
CC=[];
comparison=4;
pval = fun_compare_tastes(1,comparison, binnedC, neurons, tastes);
[CC(:,1),CC(:,2)]=find(pval<0.05);
q=1;
for i=1:size(CC,1)
    if ~isempty(intersect(CC(i,1),delay_resp_neurons))
        temp(q,1) = CC(i,1);
        temp(q,2) = CC(i,2);
        q=q+1;
    end
end
CC = temp;
plot_compare_tastes(comparison,plotbins, plotbinnedC, CC)
%%
figure
plot_taste_comparison_indiv_neuron(32,plotbins, plotbinnedC)

%% Significance between correct and error trials
pval=[];
temp=[];
CC=[];
comparison =4;
pval = fun_compare_tastes(2,comparison, binnedC3, neurons, tastes);
[CC(:,1),CC(:,2)]=find(pval<0.05);
q=1;
for i=1:size(CC,1)
    if ~isempty(intersect(CC(i,1),Bonly))
        temp(q,1) = CC(i,1);
        temp(q,2) = CC(i,2);
        q=q+1;
    end
end
CC = temp;
plot_compare_tastes(comparison,plotbins, plotbinnedC, CC)

%% plotting correct vs error
x=randi(length(CC));
figure
hold on
plot(plotbins(1:end-1), nanmean(plotbinnedC3(:,:,1,CC(x,1)),1))
plot(plotbins(1:end-1), nanmean(plotbinnedC3(:,:,3,CC(x,1)),1))
y=ylim;
line([0 0],[0 y(2)],'Color','r','LineStyle','--')
rectangle('Position',[2 0 2 y(2)],'EdgeColor', 'none', 'FaceColor', [.8, .8, .8, 0.5])
legend('Correct L','Error L')
title(append('Neuron ',num2str(CC(x,1))))
%% Plot decision Left or Right
% figure
% subplot(3,1,1)
% plot_individual_neuron(2, 48, plotbins, plotbinnedC3)
% legend('Left','Right','Location','Northwest')
% subplot(3,1,2)
% plot_individual_neuron(2, 123, plotbins, plotbinnedC3)
% subplot(3,1,3)
% plot_licking(summary)
plot_compare_L_R(plotbins, plotbinnedC3, CC(1:20))
% figure
% plot_taste_comparison_indiv_neuron(164,plotbins, plotbinnedC)
 
 %% run kruskal wallis (non-parametric ANOVA) across all tastes on different bins
 %plot significant bins/neurons on a bar graph
 %this compares response bin across tastes, does not compare to baseline
 t = categorical({'Suc','S\_75\_25','S\_55\_45','S\_45\_55','S\_25\_75','NaCl'});
%   t = categorical({'Suc','S\_55\_45','S\_45\_55','NaCl'});
comparison=3;
plot_max_response(binnedC, comparison, Aonly,t)

%% linear fit for mixture concentration
comparison=3;
% x=[100; 55; 45; 0];
rsqu=[];
x=[100; 75; 55; 45; 25; 0];

for i=1:length(taste_resp_neurons)
y=[];
md1=[];
y = squeeze(nanmean(binnedC(:,comparison,:,taste_resp_neurons(i)))); %find mean binnedC for each bin
mdl = fitlm(x,y);
rsqu(i)=mdl.Rsquared.Adjusted;

end

mixture_coding=[];
mixture_coding = find(rsqu>.7 | rsqu<-.7);
length(mixture_coding)
length(find(rsqu>.7 ))
length(find(rsqu<-.7))
%% plot all mixture coding
figure

for i=1:length(mixture_coding)
    subplot(5,4,i)
    y = squeeze(nanmean(binnedC(:,comparison,:,taste_resp_neurons(mixture_coding(i)))));
    p = polyfit(x,y,1); 
    f = polyval(p,x); 
    plot(x,y,'o',x,f,'-') 
    title(append('Neuron ',num2str(taste_resp_neurons(mixture_coding(i)))))
% legend('data','linear fit')
end
%%
figure


% for i=1:length(mixture_coding)
    subplot(2,1,1)
    y = squeeze(nanmean(binnedC(:,comparison,:,44)));
    p = polyfit(x,y,1); 
    f = polyval(p,x); 
    plot(x,y,'o',x,f,'-') 
    legend('data','linear fit')
    title(append('Neuron ',num2str(44)))
     xlabel('Sucrose Concentration')
    ylabel('Mean spiking taste response bin')
     subplot(2,1,2)
    y = squeeze(nanmean(binnedC(:,comparison,:,497)));
    p = polyfit(x,y,1); 
    f = polyval(p,x); 
    plot(x,y,'o',x,f,'-') 
    title(append('Neuron ',num2str(497)))
    xlabel('Sucrose Concentration')
    ylabel('Mean spiking taste response bin')

set(gcf, 'Renderer', 'painters');

%% mixture coding
comparison=3;
% x=[100; 55; 45; 0];
rsqu=[];
x=[100; 75; 55; 45; 25; 0];
for i=1:length(taste_resp_neurons)
y=[];
md1=[];
y = squeeze(nanmean(binnedC(:,comparison,:,taste_resp_neurons(i)))); %find mean binnedC for each bin
mdl = fitlm(x,y);
rsqu(i)=mdl.Rsquared.Ordinary;
end
%%
figure
mixture_coding = find(rsqu>.6 | rsqu<-.6);

for i=1:length(mixture_coding)
    subplot(5,6,i)
    y = squeeze(nanmean(binnedC(:,comparison,:,taste_resp_neurons(mixture_coding(i)))));
    p = polyfit(x,y,2); 
    f = polyval(p,x); 
    plot(x,y,'o',x,f,'-') 
    title(append('Neuron ',num2str(taste_resp_neurons(mixture_coding(i)))))
% legend('data','linear fit')
end

%% mixture suppression quadratic


x=[100; 75; 55; 45; 25; 0];
for i=1:length(taste_resp_neurons)
y = squeeze(nanmean(binnedC(:,comparison,:,taste_resp_neurons(i))));

p=polyfit(x,y,2); 
yfit = polyval(p,x); 
 yresid = y - yfit;
 SSresid = sum(yresid.^2);
SStotal = (length(y)-1) * var(y);
rsq(i) = 1 - SSresid/SStotal
 
% quadrsqu(i) = rsqu;
end
mixture_suppression=[]
mixture_suppression = find(rsq>.7 | rsq<-.7);
%%

figure
for i=1:length(mixture_suppression)
    subplot(6,8,i)
    y = squeeze(nanmean(binnedC(:,comparison,:,taste_resp_neurons(mixture_suppression(i)))));
    p = polyfit(x,y,2); 
    f = polyval(p,x); 
    plot(x,y,'o',x,f,'-') 
    title(append('Neuron ',num2str(taste_resp_neurons(mixture_suppression(i)))))
% legend('data','linear fit')
end
%% mixture suppression linear
comparison=3;
% x=[100; 55; 45; 0];
rsqu=[];
x=[100; 75; 55];
i=1;
for i=1:length(Aonly)
y=[];
md1=[];
y = squeeze(nanmean(binnedC(:,comparison,:,Aonly(i)))); %find mean binnedC for each bin
y2(1) = mean(y([1 6]));
y2(2) = mean(y([2 5]));
y2(3) = mean(y([3 4]));
mdl = fitlm(x,y2);
rsqu(i)=mdl.Rsquared.Adjusted;
end

mixture_supp_linear_pos = find(rsqu>.8);
mixture_supp_linear_neg = find(rsqu<-.8);
mixture_supp_linear = [mixture_supp_linear_pos mixture_supp_linear_neg];
length(mixture_supp_linear_pos)
length(mixture_supp_linear_neg)
%%
figure
for i=1:length(mixture_supp_linear)
    subplot(5,3,i)
    y = squeeze(nanmean(binnedC(:,comparison,:,Aonly(mixture_supp_linear(i)))));
    y2(1) = mean(y([1 6]));
    y2(2) = mean(y([2 5]));
    y2(3) = mean(y([3 4]));
    p = polyfit(x,y2,1); 
    f = polyval(p,x); 
    plot(x,y2,'o',x,f,'-') 
    title(append('Neuron ',num2str(Aonly(mixture_supp_linear(i)))))
% legend('data','linear fit')
end
%% plot traces
figure(124);
subplot(3,1,1)
ylimit = plot_individual_neuron(1,141, plotbins, plotbinnedC2);
hold on
rectangle('Position',[0 0 2 ylimit],'EdgeColor', 'none', 'FaceColor', [.8, .8, .8, 0.5])
hold off
subplot(3,1,2)
ylimit = plot_individual_neuron(1,733, plotbins, plotbinnedC2);
hold on
rectangle('Position',[0 0 2 ylimit],'EdgeColor', 'none', 'FaceColor', [.8, .8, .8, 0.5])
hold off
subplot(3,1,3)
plot_licking(summary)
set(gcf, 'Renderer', 'painters');
%% plot examples mixture suppression linear
figure(555)
subplot(2,1,1)
 y = squeeze(nanmean(binnedC(:,comparison,:,141)));
    y2(1) = mean(y([1 6]));
    y2(2) = mean(y([2 5]));
    y2(3) = mean(y([3 4]));
    p = polyfit(x,y2,1); 
    f = polyval(p,x); 
    plot(x,y2,'o',x,f,'-') 
    title(append('Neuron ',num2str(141)))
    xlabel('Sucrose Concentration')
    ylabel('Mean spiking taste response bin')
    legend('data','fit')
  subplot(2,1,2)  
     y = squeeze(nanmean(binnedC(:,comparison,:,733)));
    y2(1) = mean(y([1 6]));
    y2(2) = mean(y([2 5]));
    y2(3) = mean(y([3 4]));
    p = polyfit(x,y2,1); 
    f = polyval(p,x); 
    plot(x,y2,'o',x,f,'-') 
    title(append('Neuron ',num2str(733)))
   xlabel('Sucrose Concentration')
    ylabel('Mean spiking taste response bin')

set(gcf, 'Renderer', 'painters');
%%
 comparison =4;
 p=[];
 p=NaN(length(neurons),1); %initialize p value
 X = binnedC(:,comparison,:,:) ; %look at a specific time bin of binnedC
 X = squeeze(X);
 for i =1:size(X,3)
     p(i,1)= kruskalwallis(squeeze(X(:,:,i)),[],'off') ; %testing for each bin
 end
b=find(p<0.05);%find neuron and which bin they are significant for plotting
plot_max_response(binnedC, comparison, b,t)


%% plot averaged population activity

popavg1corr = nanmean(plotbinnedC4(:,:,1:2,Aonly),[1 3 4]);
popavg2corr = nanmean(plotbinnedC4(:,:,1:2,Bonly),[1 3 4]);
popavg3corr = nanmean(plotbinnedC4(:,:,1:2,Conly),[1 3 4]);
popavg1err = nanmean(plotbinnedC4(:,:,3:4,Aonly),[1 3 4]);
popavg2err = nanmean(plotbinnedC4(:,:,3:4,Bonly),[1 3 4]);
popavg3err = nanmean(plotbinnedC4(:,:,3:4,Conly),[1 3 4]);
figure
subplot(4,1,1)
hold on
plot(plotbins(1:end-1),popavg1corr)
plot(plotbins(1:end-1),popavg1err)

title('Population Averages - Taste')
y=ylim;
line([0 0],[0 y(2)],'Color','r','LineStyle','--')
rectangle('Position',[0 0 2 y(2)],'EdgeColor', 'none', 'FaceColor', [.8, .8, .8, 0.5])
legend('Correct','Error','Location','Northwest')
ylabel('Inferred Spiking (AU)')
hold off
subplot(4,1,2)
hold on
plot(plotbins(1:end-1),popavg2corr)
plot(plotbins(1:end-1),popavg2err)
title('Delay')
y=ylim;
line([0 0],[0 y(2)],'Color','r','LineStyle','--')
rectangle('Position',[2 0 2 y(2)],'EdgeColor', 'none', 'FaceColor', [.8, .8, .8, 0.5])
ylabel('Inferred Spiking (AU)')
hold off
subplot(4,1,3)
hold on
plot(plotbins(1:end-1),popavg3corr)
plot(plotbins(1:end-1),popavg3err)
title('Choice')
y=ylim;
line([0 0],[0 y(2)],'Color','r','LineStyle','--')
rectangle('Position',[4 0 2 y(2)],'EdgeColor', 'none', 'FaceColor', [.8, .8, .8, 0.5])
ylabel('Inferred Spiking (AU)')
xlabel('Time (sec)')
hold off
subplot(4,1,4)
plot_licking(summary)
% ylim([0 .012])

%  %%
%  [r,tbl,stats] = kruskalwallis(X(:,:,44))
%  intersect(a,c)
%  tbl2=ANOVA(
%  d=find(pval(6,:)<0.05);
%  setdiff(c,d)
%   setdiff(d,c)
  
 
%% spatial mapping
[correlationmatrix,center,sortedD] = Spatial_organization_CNMFE_Josh(neuron);
a=find(sortedD(:,1)==0); %remove pairs with the same neuron
sortedD(a,:)=[];
mean(sortedD(a,1))


avg_cor = mean(sortedD(:,2));
sdev = (std(sortedD(:,2)));
twosdev = (std(sortedD(:,2)))*2;
EE = avg_cor+twosdev;
FF=find(sortedD(:,2)<-EE);
GG=find(sortedD(:,2)>EE);
HH = [ GG' FF'];
HH=HH';
II=[];
II = 1:length(sortedD(:,2));
JJ = setdiff(II,HH);

histobins = -.4:.005:.4;
%%
figure
histogram(sortedD(JJ,2),histobins)
y=ylim;
hold on
histogram(sortedD(HH,2),histobins)
line([-EE -EE], [0 y(2)],'LineStyle','--','Color','r');
line([EE EE], [0 y(2)],'LineStyle','--','Color','r');
xlim([-.4 .4])
xlabel('Correlation')
ylabel('Number of Pairs')
title('Distribution of correlations between neuron pairs')
legend('Mean distance = 301.7 pixels','Mean distance = 276.3 pixels','Two standard deviations','Location','Northwest')
%%
[h,p]=ttest2(sortedD(HH,1),sortedD(JJ,1));
fprintf(append('mean correlation= ',num2str(avg_cor)))
fprintf(append('\nstdev correlation= ',num2str(sdev)))
fprintf(append('\nmean distance high correlation= ',num2str(mean(sortedD(HH,1)))))
fprintf(append('\nmean distance low correlation= ',num2str(mean(sortedD(JJ,1)))))
fprintf(append('\nP value= ',num2str(p)))
fprintf('\n')

%%
[A,loc1]=ismember(sortedD(:,3),taste_resp_neurons);
B = find(loc1);
[A,loc2]=ismember(sortedD(:,4),taste_resp_neurons);
C = find(loc2);
CC = intersect(B,C);
taste_spatial=sortedD(CC,:);
X1=abs(taste_spatial(:,2));
X2=abs(sortedD(:,2));
% histogram(taste_spatial(:,2))
hold on 
% histogram(sortedD(:,2),histobins)
X11 = mean(X1);
sem1 = std(X1)/sqrt(length(X1));
X22 = mean(X2);
sem2 = std(X2)/sqrt(length(X2));
X3 = [X11 X22];
sem3 = [sem1 sem2];
labels = categorical({'Taste responsive neurons','All neurons'});
labels = reordercats(labels,{'Taste responsive neurons','All neurons'});
figure
bar(labels, X3)
hold on
er = errorbar(X3,sem3);    
er.Color = [0 0 0];                            
er.LineStyle = 'none'; 
title('Average Correlations between pairs')
ylabel('|Average Correlation|')
[h,p]=ttest2(X1,X2)
%%
Q = neuron.C(choice_resp_neurons,:);
QQ = linkage(Q,'Ward');
% QQQ = cluster(QQ, 'Cutoff', 3, 'Depth',5);
QQQ = cluster(QQ, 'MaxClust', 6);
figure
dendrogram(QQ)
title('Hierarchical Clustering of choice responsive neurons')
neuron_number = 1:length(neurons);
clust1_neurons = choice_resp_neurons(find(QQQ == 1));
clust2_neurons = choice_resp_neurons(find(QQQ == 2));
clust3_neurons = choice_resp_neurons(find(QQQ == 3));
clust4_neurons = choice_resp_neurons(find(QQQ == 4));
clust5_neurons = choice_resp_neurons(find(QQQ == 5));
clust6_neurons = choice_resp_neurons(find(QQQ == 6));
% clust4_neurons = neuron_number(find(QQQ == 4));
% clust5_neurons = neuron_number(find(QQQ == 5));
t = categorical({'Suc','S\_75\_25','S\_55\_45','S\_45\_55','S\_25\_75','NaCl'});
%%

figure
hold on
plot(plotbins(1:end-1),nanmean(plotbinnedC2(:,:,clust1_neurons), [1 3]))
plot(plotbins(1:end-1),nanmean(plotbinnedC2(:,:,clust2_neurons), [1 3]))
plot(plotbins(1:end-1),nanmean(plotbinnedC2(:,:,clust3_neurons), [1 3]))
plot(plotbins(1:end-1),nanmean(plotbinnedC2(:,:,clust4_neurons), [1 3]))
plot(plotbins(1:end-1),nanmean(plotbinnedC2(:,:,clust5_neurons), [1 3]))
plot(plotbins(1:end-1),nanmean(plotbinnedC2(:,:,clust6_neurons), [1 3]))

y=ylim;
line([0 0], [0 y(2)],'LineStyle','--','Color','r');
rectangle('Position',[4 0 2 y(2)],'EdgeColor', 'none', 'FaceColor', [.8, .8, .8, 0.5])
legend('Cluster 1','Cluster 2','Cluster 3','Cluster 4','Cluster 5','Cluster 6','Location','Northwest')
xlabel('Time (sec)')
ylabel('Inferred Spiking (AU)')
title('Population Averages - Choice Neurons')
% plot(plotbins(1:end-1),nanmean(plotbinnedC2(:,:,clust4_neurons), [1 3]))
% plot(plotbins(1:end-1),nanmean(plotbinnedC2(:,:,clust5_neurons), [1 3]))
%%
R = neuron.C(:,:);
RR = linkage(R,'Ward');
% QQQ = cluster(QQ, 'Cutoff', 3, 'Depth',5);
RRR = cluster(RR, 'MaxClust', 6);
figure
dendrogram(RR)
title('Hierarchical Clustering of taste responsive neurons')
neuron_number = 1:length(neurons);
clust1_neurons = taste_resp_neurons(find(RRR == 1));
clust2_neurons = taste_resp_neurons(find(RRR == 2));
clust3_neurons = taste_resp_neurons(find(RRR == 3));
clust4_neurons = taste_resp_neurons(find(RRR == 4));
clust5_neurons = taste_resp_neurons(find(RRR == 5));
clust6_neurons = taste_resp_neurons(find(RRR == 6));
% clust4_neurons = neuron_number(find(QQQ == 4));
% clust5_neurons = neuron_number(find(QQQ == 5));

%%

figure
hold on
plot(plotbins(1:end-1),nanmean(plotbinnedC2(:,:,clust1_neurons), [1 3]))
plot(plotbins(1:end-1),nanmean(plotbinnedC2(:,:,clust2_neurons), [1 3]))
plot(plotbins(1:end-1),nanmean(plotbinnedC2(:,:,clust3_neurons), [1 3]))
plot(plotbins(1:end-1),nanmean(plotbinnedC2(:,:,clust4_neurons), [1 3]))
plot(plotbins(1:end-1),nanmean(plotbinnedC2(:,:,clust5_neurons), [1 3]))
plot(plotbins(1:end-1),nanmean(plotbinnedC2(:,:,clust6_neurons), [1 3]))

y=ylim;
line([0 0], [0 y(2)],'LineStyle','--','Color','r');
rectangle('Position',[4 0 2 y(2)],'EdgeColor', 'none', 'FaceColor', [.8, .8, .8, 0.5])
legend('Cluster 1','Cluster 2','Cluster 3','Cluster 4','Cluster 5','Cluster 6','Location','Northwest')
xlabel('Time (sec)')
ylabel('Inferred Spiking (AU)')
title('Population Averages - Taste Neurons')
% plot(plotbins(1:end-1),nanmean(plotbinnedC2(:,:,clust4_neurons), [1 3]))
% plot(plotbins(1:end-1),nanmean(plotbinnedC2(:,:,clust5_neurons), [1 3]))
%% cluster all neurons

R = neuron.C(:,:);
RR = linkage(R,'Ward');
% QQQ = cluster(QQ, 'Cutoff', 3, 'Depth',5);
RRR = cluster(RR, 'MaxClust', 6);
figure
dendrogram(RR)
title('Hierarchical Clustering of taste responsive neurons')
neuron_number = 1:length(neurons);
clust1_neurons = (find(RRR == 1));
clust2_neurons = (find(RRR == 2));
clust3_neurons = (find(RRR == 3));
clust4_neurons = (find(RRR == 4));
clust5_neurons = (find(RRR == 5));
clust6_neurons = (find(RRR == 6));
% clust4_neurons = neuron_number(find(QQQ == 4));
% clust5_neurons = neuron_number(find(QQQ == 5));

%%

figure
hold on
plot(plotbins(1:end-1),nanmean(plotbinnedC2(:,:,clust1_neurons), [1 3]))
plot(plotbins(1:end-1),nanmean(plotbinnedC2(:,:,clust2_neurons), [1 3]))
plot(plotbins(1:end-1),nanmean(plotbinnedC2(:,:,clust3_neurons), [1 3]))
plot(plotbins(1:end-1),nanmean(plotbinnedC2(:,:,clust4_neurons), [1 3]))
plot(plotbins(1:end-1),nanmean(plotbinnedC2(:,:,clust5_neurons), [1 3]))
plot(plotbins(1:end-1),nanmean(plotbinnedC2(:,:,clust6_neurons), [1 3]))

y=ylim;
line([0 0], [0 y(2)],'LineStyle','--','Color','r');
% rectangle('Position',[4 0 2 y(2)],'EdgeColor', 'none', 'FaceColor', [.8, .8, .8, 0.5])
legend('Cluster 1','Cluster 2','Cluster 3','Cluster 4','Cluster 5','Cluster 6','Location','Northwest')
xlabel('Time (sec)')
ylabel('Inferred Spiking (AU)')
title('Population Averages - Taste Neurons')
% plot(plotbins(1:end-1),nanmean(plotbinnedC2(:,:,clust4_neurons), [1 3]))
% plot(plotbins(1:end-1),nanmean(plotbinnedC2(:,:,clust5_neurons), [1 3]))
%%
figure
first =7:9;
second=10:12;
subplot(4,1,1)
hold on
plot(plotbins(1:end-1),nanmean(plotbinnedC4(:,:,first,clust1_neurons), [1 3 4]))
plot(plotbins(1:end-1),nanmean(plotbinnedC4(:,:,second,clust1_neurons), [1 3 4]))
% plot(plotbins(1:end-1),nanmean(plotbinnedC(:,:,3,clust1_neurons), [1 3 4]))
y=ylim;
line([0 0], [0 y(2)],'LineStyle','--','Color','r');
rectangle('Position',[4 0 2 y(2)],'EdgeColor', 'none', 'FaceColor', [.8, .8, .8, 0.5])
title('Cluster 1 Population Average')
legend('Left Error', 'Right Error','Location','Northwest')
ylabel('Inferred Spiking (AU)')

subplot(4,1,2)
hold on
plot(plotbins(1:end-1),nanmean(plotbinnedC4(:,:,first,clust2_neurons), [1 3 4]))
plot(plotbins(1:end-1),nanmean(plotbinnedC4(:,:,second,clust2_neurons), [1 3 4]))
% plot(plotbins(1:end-1),nanmean(plotbinnedC(:,:,3,clust2_neurons), [1 3 4]))
y=ylim;
line([0 0], [0 y(2)],'LineStyle','--','Color','r');
rectangle('Position',[4 0 2 y(2)],'EdgeColor', 'none', 'FaceColor', [.8, .8, .8, 0.5]);
ylabel('Inferred Spiking (AU)')
title('Cluster 2 Population Average')

subplot(4,1,3)
hold on
plot(plotbins(1:end-1),nanmean(plotbinnedC4(:,:,first,clust4_neurons), [1 3 4]))
plot(plotbins(1:end-1),nanmean(plotbinnedC4(:,:,second,clust4_neurons), [1 3 4]))
% plot(plotbins(1:end-1),nanmean(plotbinnedC(:,:,3,clust2_neurons), [1 3 4]))
y=ylim;
line([0 0], [0 y(2)],'LineStyle','--','Color','r');
rectangle('Position',[4 0 2 y(2)],'EdgeColor', 'none', 'FaceColor', [.8, .8, .8, 0.5]);
ylabel('Inferred Spiking (AU)')
title('Cluster 3 Population Average')

subplot(4,1,4)
hold on
plot(plotbins(1:end-1),nanmean(plotbinnedC4(:,:,first,clust5_neurons), [1 3 4]))
plot(plotbins(1:end-1),nanmean(plotbinnedC4(:,:,second,clust5_neurons), [1 3 4]))
% plot(plotbins(1:end-1),nanmean(plotbinnedC(:,:,3,clust3_neurons), [1 3 4]))
y=ylim;
rectangle('Position',[4 0 2 y(2)],'EdgeColor', 'none', 'FaceColor', [.8, .8, .8, 0.5]);
line([0 0], [0 y(2)],'LineStyle','--','Color','r');
ylabel('Inferred Spiking (AU)')
title('Cluster 4 Population Average')
xlabel('Time (sec)')
ylabel('Inferred Spiking (AU)')
%%
X1=abs(taste_spatial(:,1));
X2=abs(sortedD(:,1));
% histogram(taste_spatial(:,2))
hold on 
% histogram(sortedD(:,2),histobins)
X11 = mean(X1);
sem1 = std(X1)/sqrt(length(X1));
X22 = mean(X2);
sem2 = std(X2)/sqrt(length(X2));
X3 = [X11 X22];
sem3 = [sem1 sem2];
labels = categorical({'Taste responsive neurons','All neurons'});
labels = reordercats(labels,{'Taste responsive neurons','All neurons'});
figure
bar(labels, X3)
hold on
er = errorbar(X3,sem3);    
er.Color = [0 0 0];                            
er.LineStyle = 'none'; 
title('Average distance between pairs')
ylabel('Average Distance (pixels)')
[h,p]=ttest2(X1,X2)


%%
figure(123);
subplot(2,1,1)
Z= mean(plotbinnedC2(:,:,41));
ylimit = max(Z);

% rectangle('Position',[0 0 2 ylimit],'EdgeColor', 'none', 'FaceColor', [.8, .8, .8, 0.5])
bar(plotbins(1:end-1), Z);
% line([0 0],[0 ylimit],'Color','r','LineStyle','--')



hold off
subplot(2,1,2)
plot_licking(summary)

%% plot bar chart for maximum response bin

%% plot responses in peak window for individual neurons
figure
subplot(2,1,1)
plot_max_response_individual(binnedC, 3, 44,t)
subplot(2,1,2)
plot_max_response_individual(binnedC, 3, 230,t)
% subplot(3,1,3)
% plot_max_response_individual(binnedC, 5, 44,t)

%% plotting from binnedC2
figure
r1=mean(binnedC2(:,:,1,:),4);
err12 = std(r1);
r2=mean(r1);
r3=mean(binnedC2(:,:,2,:),4);
err22 = std(r2);
r4=mean(r3);
plot(bins,r2)
hold on
plot(bins,r4)
legend('Left','Right','Location','Northwest')
xlabel('Time (sec)')
ylabel('Inferred Spikes (AU)')
title('All neurons')


%% binnedC3
%dim3 1 -  correct, 2 -  error
bins = -4:.5:6; %bin size
z=1;
p=1;
binnedC3 = NaN(length(trial),length(bins),2,length(neurons));
for z = 1:length(neurons)
    x=1;
    for p=1:length(tastes)
        
        for i=1:length(neurons(z).(tastes{p}).spikes)
            for j=1:length(bins)-1
                idx= find(neurons(z).(tastes{p}).Frames{i,1} > bins(j) & neurons(z).(tastes{p}).Frames{i,1} < bins(j+1));
                r = 1+cell2mat(neurons(z).(tastes{p}).correct_choice(i,1));
                binnedC3(x,j,r,z) = sum(neurons(z).(tastes{p}).spikes{i,1}(1,idx)); 
            end  
            x=x+1;
        
        end
        x=x+1;
        
    end
end
%% plotting from binnedC3, all neurons are averaged together
%dim3 1 -  correct, 2 -  error
figure
subplot(2,1,1)
q1=nanmean(binnedC3(:,:,2,Conly),4);
err1 = std(q1);
q2=nanmean(q1);
q3=nanmean(binnedC3(:,:,1,Conly),4);
err2 = std(q3);
q4=nanmean(q3);
plot(bins(1:5),q2)
hold on
plot(bins(1:5),q4)
legend('Correct','Error','Location','Northwest')
xlabel('Time (sec)')
ylabel('Inferred Spikes (AU)')
title('All neurons')
%% dim3 , 1=
bins = -4:.5:6; %bin size
z=1;
p=1;
binnedC4 = NaN(length(trial),length(bins),4,length(neurons));
for z = 1:length(neurons)
    x=1;
    for p=1:length(tastes)
        
        for i=1:length(neurons(z).(tastes{p}).spikes)
            if cell2mat(neurons(z).(tastes{p}).correct_choice(i,1)) ==0
                r = cell2mat(neurons(z).(tastes{p}).L_R_trial(1,1));
            end
            if cell2mat(neurons(z).(tastes{p}).correct_choice(i,1)) ==1 && cell2mat(neurons(z).(tastes{p}).L_R_trial(1,1)) == 1
                r=3;
            end
            if cell2mat(neurons(z).(tastes{p}).correct_choice(i,1)) ==1 && cell2mat(neurons(z).(tastes{p}).L_R_trial(1,1)) == 2
                r=4;
            end
            for j=1:length(bins)-1
                idx= find(neurons(z).(tastes{p}).Frames{i,1} > bins(j) & neurons(z).(tastes{p}).Frames{i,1} < bins(j+1));
                binnedC4(x,j,r,z) = sum(neurons(z).(tastes{p}).spikes{i,1}(1,idx)); 
            end  
            x=x+1;
        end
        
    end
end
%% plotting from binnedC4, correct versus error trials and r vs l, all neurons are averaged together
figure
for i=1:4
    q1=nanmean(binnedC4(:,:,i,:),4);
    q2(i,:)=nanmean(q1); 
end

plot(bins(1:5),q2)
legend('Error - Left','Error - Right','Correct-Left','Correct-Right','Location','Northwest')
xlabel('Time (sec)')
ylabel('Inferred Spikes (AU)')
title('All neurons')

%% Analysis for correct vs error trials
PP = NaN( length(neurons),4);
z=1;
for j =15:17
for i=1:length(neurons)
    PP(i,z)=ranksum(binnedC3(:,j,1,i),binnedC3(:,j,2,i));
    
end
z=z+1;
end
[EE,~]=find(PP<0.05);
e = unique(EE);
sz = round(sqrt(length(e)))+1;
figure

for j=1:length(e)
    subplot(sz,sz-1,j)
    hold on
    plot(bins, nanmean(binnedC3(:,:,1,e(j)),1))
    plot(bins, nanmean(binnedC3(:,:,2,e(j)),1))
    y=ylim;
    line([0 0],[0 y(2)],'Color','r','LineStyle','--')
    title(append('Neuron ',num2str(e(j))))
%     if j==1
%         legend('S\_55\_45','S\_45\_55','Location','Northwest')
%         legend('Suc','NaCl','Location','Northwest')
%     end
    hold off
end
% suptitle('Suc vs NaCl');
% suptitle('S\_55\_45 vs S\_45\_55');

%%
test = nanmean(binnedC3(:,:,1,:),1);
size(test)
test1 = permute(test,[4,2,1,3]);
plot(bins,mean(test1))

test2 = nanmean(binnedC3(:,:,2,:),1);
size(test2)
test3 = permute(test2,[4,2,1,3]);
figure
plot(bins,mean(test1))
hold on
plot(bins,mean(test3))
[coeff1,score1,latent1] = pca(test3);
Xcentered1 = score1*coeff1';
figure
plot3(Xcentered(1,:),Xcentered(2,:),Xcentered(3,:))
hold on
plot3(Xcentered1(1,:),Xcentered1(2,:),Xcentered1(3,:))
grid on
figure
plot3(score(1,:),score(2,:),score(3,:))
hold on
plot3(score1(1,:),score1(2,:),score1(3,:))


biplot(coeff(:,1:3),'scores',score(:,1:3));
trialSUM = sum(binnedC,1);
trialMEAN = mean(binnedC,1);

%%
for j=1:length(tastes)
    for i=1:length(neurons)
        pval(j,i)=ranksum(binnedC(:,6,j,i),binnedC(:,11,j,i));
        
    end
end

find(pval1 <0.05)
 c=find(pval(1,:)<0.05);

 d=find(pval(6,:)<0.05);
 setdiff(c,d)
  setdiff(d,c)
 figure;
 for i =1:length(tastes)
    subplot(2,3,i)
    bar(bins, trialMEAN(:,:,i,330)) %plot for a given neuron
    y(i,:) = ylim; %grab ylims for each subplot
    title((tastes{i}),'Interpreter','none');
 end
 for i =1:length(tastes) % set all y limits the same
    subplot(2,3,i);
    ylim([0 max(y(:,2))]);
 end
 
 baseline = mean(binnedC(:,3:7,:,:),2);
 avgbaseline = mean(baseline,1);
 for j=1:length(tastes)
    for i=1:length(neurons)
        pvalbaseline(j,i)=ranksum(baseline(:,1,j,i),binnedC(:,8,j,i));
        
    end
 end

  c=find(pvalbaseline(1,:)<0.05);
 d=find(pval(6,:)<0.05);
 figure
 hold on
 plot(bins,nanmean(binnedC(:,:,1,30)))
  plot(bins,nanmean(binnedC(:,:,6,30)))
 
 %% plot all tastes for a given neuron
 figure 
 for i =1:length(tastes)
    subplot(2,3,i)
    bar(bins, trialMEAN(:,:,i,158)) %plot for a given neuron
    y(i,:) = ylim; %grab ylims for each subplot
    title((tastes{i}),'Interpreter','none');
 end
 for i =1:length(tastes) % set all y limits the same
    subplot(2,3,i);
    ylim([0 max(y(:,2))]);
 end
 size(binnedC);

%% average all trials and neurons and plot inferred spiking
q=mean(binnedC,3);
size(q);
q1=mean(q,4);
size(q1);
err = std(q1);
q2=mean(q1);

errorbar(bins,q2,err)
xlabel('Time (sec)')
ylabel('Inferred Spikes (AU)')
title('All neurons')
%% bin the spike data and plot for each taste separately
close;
bins = [-5:.2:10]; %bin size
neuron2plot = 52; %choose your neuron

figure
for p=1:length(tastes)
    binnedC = zeros(length(neurons(neuron2plot).(tastes{p}).spikes),length(bins));
    numberperbin = zeros(length(neurons(neuron2plot).(tastes{p}).spikes),length(bins));
for i=1:length(neurons(neuron2plot).(tastes{p}).spikes)
    for j=1:length(bins)-1
        c=1;
        for x=1:length(neurons(neuron2plot).(tastes{p}).spikes{i,1})
            if neurons(neuron2plot).(tastes{p}).Frames{i,1}(1,x) > bins(j) && neurons(neuron2plot).(tastes{p}).Frames{i,1}(1,x) < bins(j+1)
                binnedC(i,j) = (binnedC(i,j)+neurons(neuron2plot).(tastes{p}).spikes{i,1}(1,x)); %sum C value for each bin
                numberperbin(i,j) = numberperbin(i,j)+1; %keep track of number of points per bin
                
            end
        end
    end
end
y1 = ylim;
y(p) = y1(2);

sumSpikes.(tastes{p}) = sum(binnedC);
maximum(p) = max(sumSpikes.(tastes{p}));

end
ymax = max(maximum);
for p=1:length(tastes)

subplot(3,2,p)
bar( bins,sumSpikes.(tastes{p}));
ylim([0 ymax+2]);
title(append('Neuron ', num2str(neuron2plot),', taste ',(tastes{p})),'interpreter','none');
end
% figure;

% title(append('Sum of spikes for T_1 in neuron ' , num2str(neuron2plot)));
% binnedC = binnedC./numberperbin;

% for a=1:length(sumSpikes)
%     if isnan(sumSpikes(1,a))
%         sumSpikes(1,a) =0;
%     end
% end
% figure;
% plot(bins,binnedC);
% figure;
% subplot(2,1,1);


% figure
% heatmap(binnedC);
%% plot all trials for a given neuron
figure
subplot(2,1,1);


for i =1:length(trial)
%     if trial(i).TasteID == 'T_1'
         plot(trial(i).Frames_adjusted,trial(i).traces(30,:)) %number here is which cell you are looking at
%         histogram(trial(i).spikes(30,:)) %doesnt work
%     end
   hold on
end

%% plot all trials for a given neuron in binned spike data
figure
subplot(2,1,1);
bins = [-5:.05:10]; %bin size
neuron2plot = 32; %choose your neuron
p=length(tastes);
binnedC = zeros(length(neurons(neuron2plot).(tastes{p}).spikes),length(bins));
numberperbin = zeros(length(neurons(neuron2plot).(tastes{p}).spikes),length(bins));

for p=1:length(tastes)

for i=1:length(neurons(neuron2plot).(tastes{p}).spikes)
    for j=1:length(bins)-1
        c=1;
        for x=1:length(neurons(neuron2plot).(tastes{p}).spikes{i,1})
            if neurons(neuron2plot).(tastes{p}).Frames{i,1}(1,x) > bins(j) && neurons(neuron2plot).(tastes{p}).Frames{i,1}(1,x) < bins(j+1)
                binnedC(i,j) = (binnedC(i,j)+neurons(neuron2plot).(tastes{p}).spikes{i,1}(1,x)); %sum S value for each bin
                numberperbin(i,j) = numberperbin(i,j)+1; %keep track of number of points per bin
                
            end
        end
    end
end
y1 = ylim;
y(p) = y1(2);

sumSpikestotal = sum(binnedC);

% maximum(p) = max(sumSpikes.(tastes{p}));

end
bar( bins,sumSpikestotal); 
title(append('Neuron ', num2str(neuron2plot)));
ylabel('Spikes');
x=xlim;

%% plot individual spikes for all trials for a given neuron
bins = [-5:.05:10]; %bin size
neuron2plot = 20; %choose your neuron
for i=1:length(trial)
   for j=1:length(bins)-1
       for x=1:length(trial(i).imaging_frames)
       if  trial(i).imaging_frames(x) > bins(j) && trial(i).imaging_frames(x) < bins(j+1)
            if ~isequal(trial(i).spikes(neuron2plot,x),0)
                binneddata(i,x) = trial(i).imaging_frames(x);
            else
               
                binneddata(i,x) = 0;
            end
       end
       end
   end
end
%% plot licking events for all trials on same plot as spikes
neuron2plot = 3; %choose your neuron
figure;
subplot(2,1,2);
for i =1:length(binneddata(:,1))
    for j=1:length(binneddata(1,:))
        if ~isequal(binneddata(i,j),0) 
            line([binneddata(i,j) binneddata(i,j)], [i-1 i]) 
        end
    end
end
title(append('Individual spikes -','Neuron ', num2str(neuron2plot)));
ylabel('Trial #')
x=xlim;
y=ylim;

%% plot licks
subplot(2,1,2);
% figure;
for i=1:length(trial)
    j = ones(1,length(trial(i).centSp))*i;
    h1=scatter(trial(i).centSp,j,6,'filled','r'); %central licks are red
    
     hold on
end

for i=1:length(trial)
    j = ones(1,length(trial(i).LeftSp))*i;
    h2=scatter(trial(i).LeftSp,j,6,'filled','b'); %left licks are blue
     hold on
end
for i=1:length(trial)
    j = ones(1,length(trial(i).RightSp))*i;
    h3=scatter(trial(i).RightSp,j,6,'filled','g'); %right licks are green
     hold on
end
legend([h1,h2,h3],'central','left','right','Location','Northwest')
title('Licking')
ylabel('Trial #');
xlim([-4 6]);
ylim([0 length(trial)]);
%%
bins = -4:.5:6; %bin size
z=1;
p=1;
binnedC = zeros(length(neurons(z).(tastes{p}).spikes),length(bins),length(tastes),length(neurons));
for z = 1:length(neurons)
    for p=1:length(tastes)
        for i=1:length(neurons(z).(tastes{p}).spikes)
            for j=1:length(bins)-1
                idx= find(neurons(z).(tastes{p}).Frames{i,1} > bins(j) & neurons(z).(tastes{p}).Frames{i,1} < bins(j+1));
                binnedC(i,j,p,z) = sum(neurons(z).(tastes{p}).spikes{i,1}(1,idx)); %4D array, dim1=trial, dim2=bin, dim3=taste, dim4=neuron
            end  
        end
    end
end
trialSUM = sum(binnedC,1);
trialMEAN = mean(binnedC,1);
for j=1:length(tastes)
    for i=1:length(neurons)
        pval(j,i)=ranksum(binnedC(:,7,j,i),binnedC(:,15,j,i));
        
    end
end
 c=find(pval(1,:)<0.05);

 d=find(pval(6,:)<0.05);
 setdiff(c,d)
 figure;
 for i =1:length(tastes)
    subplot(2,3,i)
    bar(bins, trialMEAN(:,:,i,330)) %plot for a given neuron
    y(i,:) = ylim; %grab ylims for each subplot
    title((tastes{i}),'Interpreter','none');
 end
 for i =1:length(tastes) % set all y limits the same
    subplot(2,3,i);
    ylim([0 max(y(:,2))]);
 end
 
 baseline = mean(binnedC(:,3:7,:,:),2);
 avgbaseline = mean(baseline,1);
 for j=1:length(tastes)
    for i=1:length(neurons)
        pvalbaseline(j,i)=ranksum(baseline(:,1,j,i),binnedC(:,8,j,i));
        
    end
 end

  c=find(pvalbaseline(1,:)<0.05);

 d=find(pval(6,:)<0.05);
 figure
  for i =1:length(tastes)
    subplot(2,3,i)
    bar(bins, trialMEAN(:,:,i,158)) %plot for a given neuron
    y(i,:) = ylim; %grab ylims for each subplot
    title((tastes{i}),'Interpreter','none');
 end
 for i =1:length(tastes) % set all y limits the same
    subplot(2,3,i);
    ylim([0 max(y(:,2))]);
 end
 size(binnedC)

%% average all trials and neurons and plot inferred spiking
q=mean(binnedC,3);
size(q);
q1=mean(q,4);
size(q1);
err = std(q1);
q2=mean(q1);

errorbar(bins,q2,err)
xlabel('Time (sec)')
ylabel('Inferred Spikes (AU)')
title('All neurons')

%% distribution of neurons firing
figure;
hold on
l=[];
for i =1:2:length(bins)-1
r=mean(binnedC,3);
size(r);
r1=r(:,i,:,:);
size(r1);
err = std(r1,1);
r2=mean(r1,1);
size(r2);
r3 = reshape(r2,[],477);
spikingcells = find(r3);
% figure;
% scatter(1:length(spikingcells),r3(spikingcells))
% xlabel('Neuron')
% ylabel('Avg Inferred Spikes Preceeding lick(AU)')
% title('Distribution of Spiking')
%  subplot(6,4,i)

[h,stats] = cdfplot(r3(spikingcells));
% histogram(r3(spikingcells), 'Normalization','cdf')
ylim([.75 1])
title('CDF For all bins')

end

