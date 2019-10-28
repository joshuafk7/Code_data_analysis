% Imaging analysis for GC project
%% rigister your image
 %edit regtiff
%% load imaging data and ROI
load('D:\Behavior\Discrimination\Imaging\JK062\190925\cropped2\MC\neuron.mat')
% neuron = results;
clearvars -except neuron

%% load the event
file = dir('*.rhd');
[tastes,unit,data,trial,summary] = process_intan_v2_josh(file.name);



%% extract traces and spikes based on index of frames in C from CNMF   
for i =1:length(trial)
    trial(i).traces = neuron.C(:,trial(i).Frame_index);
    trial(i).spikes = full(neuron.S(:,trial(i).Frame_index));
end

%% reshape trials to neurons

neurons = trial2neuron5tastant_Josh(trial,tastes);

%% spatial mapping
[correlationmatrix,center,sortedD] = Spatial_organization_CNMFE_Josh(neuron);


%% bin the spike data and plot for each taste separately
bins = [-5:.5:10]; %bin size
neuron2plot = 9; %choose your neuron

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
       for x=1:length(trial(i).Frames)
       if  trial(i).Frames(x) > bins(j) && trial(i).Frames(x) < bins(j+1)
            if ~isequal(trial(i).spikes(neuron2plot,x),0)
                binneddata(i,x) = trial(i).Frames(x);
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

%plot licks
subplot(2,1,2);

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
legend([h1,h2,h3],'central','left','right')
title('Licking')
ylabel('Trial #');
xlim([-5 10]);
ylim([0 length(trial)]);

%%
% %% average resp of a given neuron across trial types
% %  for j = 1:100
% % neurons(j).centSp = [];
% %     for i =1:length(trial)
% %         for x=1:length(trial(i).centSp)
% %             
% %         neurons(j).centSp(i) =  neurons(j).centSp(i) trial(i).centSp(x);
% %         end
% %     end
% %  end
% 
% 
% 
% %% plot a given neuron traces for a given taste
% % figure
% % subplot(2,1,1);
% neurontoplot = 26; %choose which neuron to plot
% % for i = 1:length(neurons(neurontoplot).T_1) %loop over trials for T_1%plotting T_1 (ex sucrose)
% %     for x = 1:length(neurons(neurontoplot).T_1{i,1}) %loop over C values
% %         plot(neurons(neurontoplot).T_1{i,2},neurons(neurontoplot).T_1{i,1}) %%x here refers to frames for x axis corresponding to C values
% %         hold on
% %     end
% % end
% % title(append('All traces for T_1 in neuron ' , num2str(neurontoplot)));
% % hold off
% % figure
% % for i = 1:length(neurons(26).T_7) %plotting T_7 (ex NaCl)
% %     for x = 1:length(neurons(26).T_7{i,1})
% %         plot(neurons(26).T_7{i,2},neurons(26).T_7{i,1})
% %         hold on
% %     end
% % end
% 
% %mean response to a given taste for a given neuron
% % binning the values of C and then averaging when a given bin has multiple
% % C values within it
% 
% % not sure if this is correct, but plot of avgC looks like an average of
% % all of the values from plot in previous section
% 
% bins = [-5:.25:10]; %bin size
% binnedC = zeros(length(neurons(neurontoplot).T_1),length(bins));
% numberperbin = zeros(length(neurons(neurontoplot).T_1),length(bins));
% 
% for i=1:length(neurons(neurontoplot).T_1)
%     for j=1:length(bins)-1
%         c=1;
%         for x=1:length(neurons(neurontoplot).T_1{i,1})
%             if neurons(neurontoplot).T_1{i,2}(1,x) > bins(j) && neurons(neurontoplot).T_1{i,2}(1,x) < bins(j+1)
%                 binnedC(i,j) = (binnedC(i,j)+neurons(neurontoplot).T_1{i,1}(1,x)); %sum C value for each bin
%                 numberperbin(i,j) = numberperbin(i,j)+1; %keep track of number of points per bin
%                 
%             end
%         end
%     end
% end
% % binnedC = binnedC./numberperbin;
% sumSpikes = sum(binnedC);
% % for a=1:length(sumSpikes)
% %     if isnan(sumSpikes(1,a))
% %         sumSpikes(1,a) =0;
% %     end
% % end
% % figure;
% % plot(bins,binnedC);
% % figure;
% subplot(2,1,1);
% 
% bar( sumSpikes);
% title(append('Sum of spikes for T_1 in neuron ' , num2str(neurontoplot)));
% % figure
% % heatmap(binnedC);
% 
% %% test each bin for significance against baseline
% for i =1:length(binnedC(:,1))
%     baseline(i,1) = sum(binnedC(i,10:20))/10;
% end
% x=1;
% for j = 41:length(binnedC(1,:))
%     p(1,x)  = signrank(baseline(:,1), binnedC(:,j));
%     x=x+1;
% end
% 
% %% generate a trial structure; I know each trial contains 60 frames; and licking data are aligned to the tone;
% dF = reshape(full(neuron.C), size(neuron.C,2),60,[]);
% for i = 1:length(trial)
%     trial(i).licks = info(i).lick - info(i).tone;
%     trial(i).trace = dF(:,:,i);
% end
% save('data.mat','Y_r','trial')
% %% get the timestamps of each frame after 5 time averaging
% for i = 1:length(trial)
%     idx = 3:5:300;
%     trial(i).framT =trial(i).Frame(idx);
% end
% %% Smooth the time-series with gaussian kerner
% for i = 1:length(trial)
%     for j = 1:size(trial(1).trace,1)
%         trial(i).traceSmooth(j,:) = gaussmooth(trial(i).trace(j,:),5,1);
%     end  
% end
% %% calculate the dF/F0; F0 is 1 s before the tone
% for i = 1:length(trial)
%     idx = find(trial(1).framT>-1 & trial(1).framT<0);
%     baseline =mean(trial(i).traceSmooth(:,idx),2);
%     trial(i).traceSmooth_dF = (trial(i).traceSmooth-repmat(baseline,1,size(trial(i).traceSmooth,2)))./repmat(baseline,1,size(trial(i).traceSmooth,2));
% end
% %% Plot Tone response
% for j =101:117
% n = j;
% figure;
% for i = 1: length(trial)
%     plot(trial(i).framT, trial(i).traceSmooth_dF(n,:),'Color',[0,0,0,0.5])
%     hold on
% end
% title(['Neuron # ', num2str(n)])
% xlim([-1.5,7])
% 
% end
% %% Plot the Sucrose trial
% % for j = 1:size(trial(1).trace,1)
%     n = 105;
%     figure;
%     for i = 1: length(trial)
%         if ~isnan(trial(i).S)
%         plot(trial(i).framT, trial(i).traceSmooth_dF(n,:),'Color',[0,0,0,0.5])
%         hold on
%         end
%     end
%     title(['Neuron # ', num2str(n)])
%     xlim([-1.5,7])
% % end
% %% Plot the Maltose trial
% % for j = 1:size(trial(1).trace,1)
%     n = 105;
%     figure;
%     for i = 1: length(trial)
%         if ~isnan(trial(i).N)
%         plot(trial(i).framT, trial(i).traceSmooth_dF(n,:),'Color',[0,0,0,0.5])
%         hold on
%         end
%     end
%     title(['Neuron # ', num2str(n)])
%     xlim([-1.5,7])
% % end
% %% Plot the quinine trial
%     n = 105;
%     figure;
%     for i = 1: length(trial)
%         if ~isnan(trial(i).CA)
%         plot(trial(i).framT, trial(i).traceSmooth_dF(n,:),'Color',[0,0,0,0.5])
%         hold on
%         end
%     end
%     title(['Neuron # ', num2str(n)])
%     xlim([-1.5,7])
%     
% %% Plot the Cyclohexamide trial
%     n = 105;
%     figure;
%     for i = 1: length(trial)
%         if ~isnan(trial(i).Q)
%         plot(trial(i).framT, trial(i).traceSmooth_dF(n,:),'Color',[0,0,0,0.5])
%         hold on
%         end
%     end
%     title(['Neuron # ', num2str(n)])
%     xlim([-1.5,7])
%     
% %% Plot the water trial
%     n = 105;
%     figure;
%     for i = 1: length(trial)
%         if ~isnan(trial(i).W)
%         plot(trial(i).framT, trial(i).traceSmooth_dF(n,:),'Color',[0,0,0,0.5])
%         hold on
%         end
%     end
%     title(['Neuron # ', num2str(n)])
%     xlim([-1.5,7])
% % end
% %% Visualize the averaged response
% % for n = 1:27
% % n =12;
% n =105;
% it = 1;
% for i = 1: length(trial)
%     if ~isnan(trial(i).S)
%         S_trace_dF(it,:) = trial(i).traceSmooth_dF(n,:);
%         it =1+it;
%     end
% end
% 
% % n =12;
% it = 1;
% for i = 1: length(trial)
%     if ~isnan(trial(i).N)
%         M_trace_dF(it,:) = trial(i).traceSmooth_dF(n,:);
%         it =1+it;
%     end
% end
% 
% it = 1;
% for i = 1: length(trial)
%     if ~isnan(trial(i).CA)
%         Q_trace_dF(it,:) = trial(i).traceSmooth_dF(n,:);
%         it =1+it;
%     end
% end
% 
% it = 1;
% for i = 1: length(trial)
%     if ~isnan(trial(i).Q)
%         Cy_trace_dF(it,:) = trial(i).traceSmooth_dF(n,:);
%         it =1+it;
%     end
% end
% 
% figure;
% plot(trial(1).framT,mean(S_trace_dF,1))
% hold on
% plot(trial(1).framT,mean(M_trace_dF,1))
% plot(trial(1).framT,mean(Q_trace_dF,1))
% plot(trial(1).framT,mean(Cy_trace_dF,1))
% xlim([-1.5,7])
% legend({'S','M','Q','Cy'})
% % end
% %% align to the tastant
% for i = 1:length(trial)
%     if ~isnan(trial(i).S)
%         trial(i).Time_Taste = trial(i).framT-trial(i).S(1);
%     elseif ~isnan(trial(i).N)
%         trial(i).Time_Taste = trial(i).framT-trial(i).N(1);
%     elseif ~isnan(trial(i).CA)
%         trial(i).Time_Taste = trial(i).framT-trial(i).CA(1);
%     elseif ~isnan(trial(i).Q)
%         trial(i).Time_Taste = trial(i).framT-trial(i).Q(1);
%     elseif ~isnan(trial(i).W)
%         trial(i).Time_Taste = trial(i).framT-trial(i).W(1);
%     end
% end
% for i = 1:length(trial)
%     idx = find(trial(i).Time_Taste>-4 & trial(i).Time_Taste<4);  % should be 50
%     if length(idx)==50
%         trial(i).T = trial(i).Time_Taste(idx);
%         trial(i).Taste = trial(i).traceSmooth_dF(:,idx);
%     elseif length(idx) ==49
%         idx(end+1) = idx(end)+1;
%         trial(i).T = trial(i).Time_Taste(idx);
%         trial(i).Taste = trial(i).traceSmooth_dF(:,idx);
%     end 
% end
% %%
% figure;
% for i = 1:length(trial)
%     T(i,:) = trial(i).T;
% end
% T = mean(T,1);
% for i = 1:length(trial)
%     trial(i).Tpro = T; % creat a proximate time for all tastant, as tastant may jitter a little bit.
% end
% %% plot the average response aligned to tastant
% n =105;
% it = 1;
% for i = 1: length(trial)
%     if ~isnan(trial(i).S)
%         S_Taste_dF(it,:) = trial(i).Taste(n,:);
%         it =1+it;
%     end
% end
% 
% % n =12;
% it = 1;
% for i = 1: length(trial)
%     if ~isnan(trial(i).N)
%         M_Taste_dF(it,:) = trial(i).Taste(n,:);
%         it =1+it;
%     end
% end
% 
% it = 1;
% for i = 1: length(trial)
%     if ~isnan(trial(i).W)
%         W_Taste_dF(it,:) = trial(i).Taste(n,:);
%         it =1+it;
%     end
% end
% 
% 
% figure;
% plot(T,mean(S_Taste_dF,1))
% hold on
% plot(T,mean(M_Taste_dF,1))
% plot(T,mean(W_Taste_dF,1))
% ylabel('dF/F')
% xlabel('Time (s)')
% xlim([-4,4])
% legend('Sucrose','Maltose','Water')
% % end
% %Plot the Sucrose trial
% % for j = 1:size(trial(1).trace,1)
% % n = 12;
% figure;
% for i = 1: size(S_Taste_dF,1)  
%     plot(T, S_Taste_dF(i,:),'Color',[0,0,0,0.2])
%     hold on
% end
% plot(T,mean(S_Taste_dF,1),'k')
% title(['Neuron # ', num2str(n), ' Sucrose'])
% xlim([-4,4])
% ylabel('dF/F')
% xlabel('Time (s)')   
% % Plot the Maltose trial    
% % n = 12;
% figure;
% for i = 1: size(M_Taste_dF,1)  
%     plot(T, M_Taste_dF(i,:),'Color',[0,0,0,0.2])
%     hold on
% end
% plot(T,mean(M_Taste_dF,1),'k')
% title(['Neuron # ', num2str(n), ' Maltose'])
% xlim([-4,4])
% ylabel('dF/F')
% xlabel('Time (s)')
% 
% figure;
% for i = 1: size(W_Taste_dF,1)  
%     plot(T, W_Taste_dF(i,:),'Color',[0,0,0,0.2])
%     hold on
% end
% plot(T,mean(W_Taste_dF,1),'k')
% title(['Neuron # ', num2str(n), ' Water'])
% xlim([-4,4])
% ylabel('dF/F')
% xlabel('Time (s)')
% %%
% save('data.mat','Y_r','trial')
% %% reorganize the data
% neuron = trial2neuron5tastant(trial);
% %% stats for each tastant
% %% statistical test
% for j = 1:length(neuron)
%     % j = 15;
%     idx = find(neuron(j).T>-1 & neuron(j).T <0);
%     T_idx1 = find(neuron(j).T>0 & neuron(j).T <3);
% %     T_idx2 = find(neuron(j).T>1 & neuron(j).T <2);
% %     T_idx3 = find(neuron(j).T>2 & neuron(j).T <3);
%     S_baseline    = mean(neuron(j).S_Taste_dF(:,idx),2);
%     S_Taste_1     = mean(neuron(j).S_Taste_dF(:,T_idx1),2);
% %     S_Taste_2   = mean(neuron(j).S_Taste_dF(:,T_idx2),2);
% %     S_Taste_3    = mean(neuron(j).S_Taste_dF(:,T_idx3),2);
%     
%     [p(1),h(1)] = ranksum(S_baseline,S_Taste_1);
%     if mean(S_Taste_1)< mean(S_baseline);
%         h(1) = 0;
%     end
% %     [p(2),h(2)] = ranksum(S_baseline,S_Taste_2);
% %     if mean(S_Taste_2)< mean(S_baseline);
% %         h(2) = 0;
% %     end
% %     [p(3),h(3)] = ranksum(S_baseline,S_Taste_3);
% %     if mean(S_Taste_3)< mean(S_baseline);
% %         h(3) = 0;
% %     end
%     M_baseline    = mean(neuron(j).M_Taste_dF(:,idx),2); % 2nd taste
%     M_Taste_1     = mean(neuron(j).M_Taste_dF(:,T_idx1),2);
%     [p(2),h(2)] = ranksum(M_baseline,M_Taste_1);
%     if mean(M_Taste_1)< mean(M_baseline);
%         h(2) = 0;
%     end
%     
%     CA_baseline    = mean(neuron(j).CA_Taste_dF(:,idx),2); % 3rd taste
%     CA_Taste_1     = mean(neuron(j).CA_Taste_dF(:,T_idx1),2);
%     [p(3),h(3)] = ranksum(CA_baseline,CA_Taste_1);
%     if mean(CA_Taste_1)< mean(CA_baseline);
%         h(3) = 0;
%     end
%     
%     Q_baseline    = mean(neuron(j).Q_Taste_dF(:,idx),2); % 4th taste
%     Q_Taste_1     = mean(neuron(j).Q_Taste_dF(:,T_idx1),2);
%     [p(4),h(4)] = ranksum(Q_baseline,Q_Taste_1);
%     if mean(Q_Taste_1)< mean(Q_baseline);
%         h(4) = 0;
%     end   
%     
%     W_baseline    = mean(neuron(j).W_Taste_dF(:,idx),2); % 4th taste
%     W_Taste_1     = mean(neuron(j).W_Taste_dF(:,T_idx1),2);
%     [p(5),h(5)] = ranksum(W_baseline,W_Taste_1);
%     if mean(W_Taste_1)< mean(W_baseline);
%         h(5) = 0;
%     end   
%     
%     
%     neuron(j).Sres = h(1);
%     neuron(j).Mres = h(2);
%     neuron(j).CAres = h(3);
%     neuron(j).Qres = h(4);
%     neuron(j).Wres = h(5);
% end
% %%
% plot_dF(116,neuron)
% 
% %% statistical test here the baseline is the 1 s before the cue; all tastant are tested together
% % for j = 1:length(neuron)
% %     % j = 15;
% %     idx = find(trial(1).framT>-1 & trial(1).framT <0);
% %     T_idx1 = find(trial(1).T>0 & trial(1).T <1);
% %     T_idx2 = find(trial(1).T>1 & trial(1).T <2);
% %     T_idx3 = find(trial(1).T>2 & trial(1).T <3);
% %     for i = 1:length(trial)
% %         baseline(i) = mean(trial(i).traceSmooth_dF(j,idx),2);
% %         Taste_1(i)    = mean(trial(i).Taste(j,T_idx1),2);
% %         Taste_2(i)    = mean(trial(i).Taste(j,T_idx2),2);
% %         Taste_3(i)    = mean(trial(i).Taste(j,T_idx3),2);
% %     end
% %     [p(1),h(1)] = ranksum(baseline,Taste_1);
% %     if mean(Taste_1)< mean(baseline);
% %         h(1) = 0;
% %     end
% %     [p(2),h(2)] = ranksum(baseline,Taste_2);
% %     if mean(Taste_2)< mean(baseline);
% %         h(2) = 0;
% %     end
% %     [p(3),h(3)] = ranksum(baseline,Taste_3);
% %     if mean(Taste_3)< mean(baseline);
% %         h(3) = 0;
% %     end
% %     neuron(j).TasteRes.p = p;
% %     neuron(j).TasteRes.h = sum(h);
% %     neuron(j).TasteResponse = sum(h);
% % end
% %%
% save('data.mat','trial','Y_r','neuron')