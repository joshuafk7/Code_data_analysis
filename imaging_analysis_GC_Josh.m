% Imaging analysis for GC project
%% rigister your image
 %edit regtiff
%% load imaging data and ROI
load('D:\Behavior\Discrimination\Imaging\JK105\200126\neuron.mat')
% neuron = results;

neuron = neuron_new2;
clearvars -except neuron
%% load the event
file = dir('*.rhd');
a = {'Suc','S_75_25','S_55_45','S_45_55','S_25_75', 'NaCl'};
% a = {'Suc','S_85_15','S_65_35','S_35_65','S_15_85', 'NaCl'};
b = [1 1 1 2 2 2];
[tastes,~,~,trial,summary] = process_intan_v4_behavior_only(file.name, a, b);
clear file a b;


%% extract traces and spikes to trial structurebased on index of frames in C from CNMF   
for i =1:length(trial)
    trial(i).traces = neuron.C(:,trial(i).Imaging_Frame_index);
    trial(i).spikes = full(neuron.S(:,trial(i).Imaging_Frame_index));
end

%% reshape trials to neurons

neurons = trial2neuron5tastant_Josh(trial,tastes);

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
plot(plotbins(1:end-1),nanmean(plotbinnedC2(:,:,clust4_neurons), [1 3]))
plot(plotbins(1:end-1),nanmean(plotbinnedC2(:,:,clust5_neurons), [1 3]))

y=ylim;
line([0 0], [0 y(2)],'LineStyle','--','Color','r');
rectangle('Position',[4 0 2 y(2)],'EdgeColor', 'none', 'FaceColor', [.8, .8, .8, 0.5])
legend('Cluster 1','Cluster 2','Cluster 4','Cluster 5','Location','Northwest')
xlabel('Time (sec)')
ylabel('Inferred Spiking (AU)')
title('Population Averages - Choice Neurons')
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
%% binned data - create 4d matrix 
%4D array called binnedC, dim1=trial, dim2=bin, dim3=taste, dim4=neuron
bins = -4:2:6; %bin size for comparison 2 second
binnedC = fun_binning(1,trial, bins, tastes, neurons);
plotbins = -4:.5:6; %smaller bin sizes for plotting
plotbinnedC = fun_binning(1,trial, plotbins, tastes, neurons);
%% binnedC2 - trial x bin x neuron
binnedC2 = fun_binning(2, trial, bins, tastes, neurons);
plotbinnedC2 = fun_binning(2,trial, plotbins, tastes, neurons);
%% binnedC3 trial x bin x correct vs error x neuron
binnedC3 = fun_binning(3, trial, bins, tastes, neurons);
plotbinnedC3 = fun_binning(3,trial, plotbins, tastes, neurons);
%% binnedC4 trial x bin x taste, correct vs error x neuron
binnedC4 = fun_binning(4, trial, bins, tastes, neurons);
plotbinnedC4 = fun_binning(4,trial, plotbins, tastes, neurons);
%% check multiple time bins against baseline for taste responsiveness
baseline =2;
comparison = 3;
taste_resp_neurons = fun_taste_responses(neurons, tastes, binnedC, binnedC2, baseline, comparison);
plot_sig_neurons(taste_resp_neurons, plotbins, plotbinnedC2,summary)

comparison = 4;
delay_resp_neurons = fun_taste_responses(neurons, tastes, binnedC, binnedC2, baseline, comparison);
plot_sig_neurons(taste_resp_neurons, plotbins, plotbinnedC2,summary)

comparison = 5;
choice_resp_neurons = fun_taste_responses(neurons, tastes, binnedC, binnedC2, baseline, comparison);
plot_sig_neurons(choice_resp_neurons, plotbins, plotbinnedC2,summary)


%% 
figure(123);
subplot(4,1,1)
ylimit = plot_individual_neuron(1,14, plotbins, plotbinnedC2);
hold on
rectangle('Position',[0 0 2 ylimit],'EdgeColor', 'none', 'FaceColor', [.8, .8, .8, 0.5])
hold off
subplot(4,1,2)
ylimit = plot_individual_neuron(1,114, plotbins, plotbinnedC2);
hold on
rectangle('Position',[2 0 2 ylimit],'EdgeColor', 'none', 'FaceColor', [.8, .8, .8, 0.5])
hold off
subplot(4,1,3)
ylimit =plot_individual_neuron(1,276, plotbins, plotbinnedC2);
hold on
rectangle('Position',[4 0 2 ylimit],'EdgeColor', 'none', 'FaceColor', [.8, .8, .8, 0.5])
hold off
subplot(4,1,4)
plot_licking(summary)


%% plot bar chart for maximum response bin
t = categorical({'Suc','S\_75\_25','S\_55\_45','S\_45\_55','S\_25\_75','NaCl'});
comparison=3;
plot_max_response(binnedC, comparison, taste_resp_neurons,t)
%% plot responses in peak window for individual neurons
figure
subplot(3,1,1)
plot_max_response_individual(binnedC, 3, 120,t)
subplot(3,1,2)
plot_max_response_individual(binnedC, 4, 114,t)
subplot(3,1,3)
plot_max_response_individual(binnedC, 5, 276,t)
 %% intersection analysis for venn diagram
 A= taste_resp_neurons;
 B=delay_resp_neurons;
 C=choice_resp_neurons;
 AB = intersect(A,B);
 AC = intersect(A,C);
 BC = intersect(B,C);
 ABC = intersect(AB,C);
 venn(1) = length(AB);
 venn(2) = length(AC);
 venn(3) = length(BC);
 venn(4) = length(ABC);


%% comparison between tastes
pval=[];
temp=[];
CC=[];
pval = fun_compare_tastes(1,comparison, binnedC, neurons, tastes);
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
plot_compare_tastes(plotbins, plotbinnedC, CC)
figure
plot_taste_comparison_indiv_neuron(164,plotbins, plotbinnedC)

%%
pval=[];
temp=[];
CC=[];
comparison =5;
pval = fun_compare_tastes(2,comparison, binnedC3, neurons, tastes);
[CC(:,1),CC(:,2)]=find(pval<0.05);
q=1;
for i=1:size(CC,1)
    if ~isempty(intersect(CC(i,1),choice_resp_neurons))
        temp(q,1) = CC(i,1);
        temp(q,2) = CC(i,2);
        q=q+1;
    end
end
CC = temp;
%%
figure
subplot(3,1,1)
plot_individual_neuron(2, 48, plotbins, plotbinnedC3)
legend('Left','Right','Location','Northwest')
subplot(3,1,2)
plot_individual_neuron(2, 123, plotbins, plotbinnedC3)
subplot(3,1,3)
plot_licking(summary)
plot_compare_L_R(plotbins, plotbinnedC3, CC)
% figure
% plot_taste_comparison_indiv_neuron(164,plotbins, plotbinnedC)
 
 %% run kruskal wallis (non-parametric ANOVA) across all tastes on different bins
 %plot significant bins/neurons on a bar graph
 %this compares response bin across tastes, does not compare to baseline
 comparison =3;
 p=NaN(length(neurons)); %initialize p value
 X = binnedC(:,comparison,:,:) ; %look at a specific time bin of binnedC
 X = squeeze(X);
 for i =1:size(X,3)
     p(i)= kruskalwallis(X(:,:,i),[],'off') ; %testing for each bin
 end
b=find(p<0.05);%find neuron and which bin they are significant for plotting
plot_max_response(binnedC, comparison, b,t)


%% plot averaged population activity

popavg1corr = nanmean(plotbinnedC4(:,:,1:6,taste_resp_neurons),[1 3 4]);
popavg2corr = nanmean(plotbinnedC4(:,:,1:6,delay_resp_neurons),[1 3 4]);
popavg3corr = nanmean(plotbinnedC4(:,:,1:6,choice_resp_neurons),[1 3 4]);
popavg1err = nanmean(plotbinnedC4(:,:,7:12,taste_resp_neurons),[1 3 4]);
popavg2err = nanmean(plotbinnedC4(:,:,7:12,delay_resp_neurons),[1 3 4]);
popavg3err = nanmean(plotbinnedC4(:,:,7:12,choice_resp_neurons),[1 3 4]);
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

 %%
 [r,tbl,stats] = kruskalwallis(X(:,:,44))
 intersect(a,c)
 
 d=find(pval(6,:)<0.05);
 setdiff(c,d)
  setdiff(d,c)
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
q1=nanmean(binnedC3(:,:,2,:),4);
err1 = std(q1);
q2=nanmean(q1);
q3=nanmean(binnedC3(:,:,1,:),4);
err2 = std(q3);
q4=nanmean(q3);
plot(bins,q2)
hold on
plot(bins,q4)
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

plot(bins,q2)
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

