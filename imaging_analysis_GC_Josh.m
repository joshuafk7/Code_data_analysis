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
b = [1 1 1 2 2 2];
[tastes,unit,data,trial,summary] = process_intan_v4_behavior_only(file.name, a, b);



%% extract traces and spikes to trial structurebased on index of frames in C from CNMF   
for i =1:length(trial)
    trial(i).traces = neuron.C(:,trial(i).Imaging_Frame_index);
    trial(i).spikes = full(neuron.S(:,trial(i).Imaging_Frame_index));
end

%% reshape trials to neurons

neurons = trial2neuron5tastant_Josh(trial,tastes);

%% spatial mapping
[correlationmatrix,center,sortedD] = Spatial_organization_CNMFE_Josh(neuron);
h= heatmap(correlationmatrix);
h.Colormap = parula;
%% binned data - create 4d matrix 
%4D array called binnedC, dim1=trial, dim2=bin, dim3=taste, dim4=neuron
bins = -4:.5:6; %bin size
z=1;
p=1;
binnedC = NaN(length(trial),length(bins),length(tastes),length(neurons));
for z = 1:length(neurons)
    x=1;
    for p=1:length(tastes)
        for i=1:length(neurons(z).(tastes{p}).spikes)
            for j=1:length(bins)-1
                idx= find(neurons(z).(tastes{p}).Frames{i,1} > bins(j) & neurons(z).(tastes{p}).Frames{i,1} < bins(j+1));
                binnedC(x,j,p,z) = sum(neurons(z).(tastes{p}).spikes{i,1}(1,idx)); 
            end  
            x=x+1;
        end
    end
end
%% binnedC2 - trial x bin x neuron
bins = -4:.5:6; %bin size
z=1;
p=1;
binnedC2 = NaN(length(trial),length(bins)-1,length(neurons));
for z = 1:length(neurons)
    x=1;
    for p=1:length(tastes)
        for i=1:length(neurons(z).(tastes{p}).spikes)
            for j=1:length(bins)-1
                idx= find(neurons(z).(tastes{p}).Frames{i,1} > bins(j) & neurons(z).(tastes{p}).Frames{i,1} < bins(j+1));
%                 r = cell2mat(neurons(z).(tastes{p}).L_R_trial(1,1));
                binnedC2(x,j,z) = sum(neurons(z).(tastes{p}).spikes{i,1}(1,idx)); 
            end 
            x=x+1;
        end
    end
end
%%
pval3=NaN(6,length(neurons),4);
x=[];
idx=1;
for z = 9:12
for j =1:6
    for i=1:length(neurons)
        pval3(j,i,idx)=ranksum(binnedC(:,7,j,i),binnedC(:,z,j,i));
        
    end
end
idx=idx+1;
end
pval3=pval3;
x=[];
for j=1:4
[row,col]=find(pval3(:,:,j)<0.05);
if j==1
    r1=row';
    c1 = col';
else
r1 = [r1 row'];
c1 = [c1 col'];
end
end
x(:,1) = r1';
x(:,2) = c1';
z = unique(x(:,2));
for j = 1:6
for i = 1:length(z)
    PP(i,j)=nanmean(binnedC(:,11,j,z(i)));
end
end
heatmap(PP)
size(intersect(a,z))
%%
%taste responsive neurons across all taste conditions
    for i=1:length(neurons)
        pval1(i)=ranksum(binnedC2(:,7,i),binnedC2(:,11,i));
        
    end
     for i=1:length(neurons)
        pval2(i)=ranksum(binnedC(:,7,1,i),binnedC(:,11,1,i));
        
    end
a=find(pval1 <0.05);
larger = find(nanmean(binnedC2(:,7,a))<nanmean(binnedC2(:,11,a)));
a = a(larger); %only choose excitatory responses
% b=find(pval2 <0.05);
for j=1:length(a)
    subplot(5,6,j) %plot all of the taste responsive neurons
    plot(bins(1:end-1),nanmean(binnedC2(:,:,a(j))))
    y=ylim;
    line([0 0],[0 y(2)],'Color','r','LineStyle','--')
    title(append('Neuron ',num2str(a(j))))
end
%% plot individual tastes for the comparison bin regardless of differences between tastes
figure
 t=categorical(tastes);
 for j=1:length(a)
     subplot(5,6,j)
    for i=1:6
        hold on
   
        bar1 = bar(t(i), nanmean(X(:,i,a(j))));
        
        errorbar(t(i),nanmean(X(:,i,a(j))),nanstd(X(:,i,a(j)))/sqrt(length(X(:,i,a(j)))));
        title(append('Neuron ',num2str(a(j))))
    end
    hold off
   
 end
%% comparison between tastes

for i=1:length(neurons)
    pval(i)=ranksum(binnedC(:,11,3,i),binnedC(:,11,4,i));
    
end
c=find(pval<0.05);
figure

for j=1:length(c)
    subplot(3,3,j)
    hold on
    plot(bins, nanmean(binnedC(:,:,1,c(j)),1))
    plot(bins, nanmean(binnedC(:,:,6,c(j)),1))
    y=ylim;
    line([0 0],[0 y(2)],'Color','r','LineStyle','--')
    title(append('Neuron ',num2str(c(j))))
    if j==1
        legend('Suc','NaCl','Location','Northwest')
    end
    hold off
end
%  figure
%  hold on
%   plot(bins, nanmean(binnedC(:,:,1,363),1))
%  plot(bins, nanmean(binnedC(:,:,6,363),1))
%  
 %%
X = binnedC(:,11,:,:) ;
X = squeeze(X);
 for i =1:size(X,3)
    p(i)= kruskalwallis(X(:,:,i),[],'off') ;
 end
b =  find(p<0.05);

 figure
 t=categorical(tastes);
 for j=1:length(b)
     subplot(5,3,j)
    for i=1:6
        hold on
   
        bar1 = bar(t(i), nanmean(X(:,i,b(j))));
        
        errorbar(t(i),nanmean(X(:,i,b(j))),nanstd(X(:,i,b(j)))/sqrt(length(X(:,i,b(j)))));
        title(append('Neuron ',num2str(b(j))))
    end
    hold off
   
 end
%  errorbar(nanmean(X(:,i,b(j))),nanstd(X(:,i,b(j)))/sqrt(length(X(:,i,b(j)))));

%% plot averaged population activity
temp=[];
temp = nanmean(binnedC2,3);
temp = nanmean(temp,1);
figure
plot(bins,temp)
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
        x=x+1;
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

