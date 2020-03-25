a = dir('*.csv');
coords=readtable(a(1).name);
%euclidean distance for jaw movement between t and t-1
eucdist = sqrt((coords.jaw_x(2:end) - coords.jaw_x(1:end-1)).^2 + (coords.jaw_y(2:end) - coords.jaw_y(end-1)).^2);
eucdist = eucdist';
c = {'Suc','S_85_15','S_65_35','S_35_65','S_15_85','NaCl'};
d=[1 1 1 2 2 2];
file = dir('*.rhd');
[~,~,~,~, summary] = process_intan_v4_behavior_only(file.name,c,d);
%% smooth data
windowSize = 5; 
b = (1/windowSize)*ones(1,windowSize);
a = 1;
distsmooth = filter(b,a,eucdist);

hold on;
plot(eucdist(1:5000));
plot(distsmooth(1:5000));
%% align jaw movments with trial structure
for i =1:length(summary.trial)
    if ~isempty(summary.trial(i).Beh_Frame_index)
    summary.trial(i).jaw_movement = eucdist(summary.trial(i).Beh_Frame_index(1,:));
    end
end
figure
plot(summary.trial(50).beh_frames,summary.trial(50).jaw_movement)


%% jaw movement no binning
p=1;
for i =1:length(trial)
     if length(trial(i).Frame_index)>495
   e(p,:) = trial(i).jaw_movement(1,1:497);
   p=p+1;
   
     end
end
q=mean(e);
%% binned jaw movement
p=1;
bins = -4:.2:6;

for i=1:length(summary.trial)
    if length(summary.trial(i).Beh_Frame_index)>496
        for j =1:length(bins)-1
            
            binned_data(p,j) = mean(summary.trial(i).jaw_movement(find(summary.trial(i).beh_frames > bins(j) & summary.trial(i).beh_frames < bins(j+1))));
            
        end
    p=p+1;    
    end
    
end
z(1,:) = mean(binned_data(1:60,:)); %mean distance for each bin across all trials
z(2,:) = mean(binned_data(61:end,:));
z(3,:) = mean(binned_data(:,:));
figure
subplot(3,1,1)
plot(bins(1:end-1),z(1,:))
hold on
plot(bins(1:end-1),z(2,:))
plot(bins(1:end-1),z(3,:))
hold off
title(append(summary.mouseID,' - Jaw Movement Averaged Across Trials'))
ylabel('Distance between t and t-1')
xlabel('Trial time aligned to first lick (sec)')
legend('first half','second half','all trials')
%% plot subset of trials on same graph
figure


for i=10:15 %which trials you want
    
      if ~isempty(trial(i).Frame_index)
    plot(trial(i).beh_frames,trial(i).jaw_movement)
    hold on
    
      end
end
%% plot jaw movement along with lick raster from lick detection
subplot(2,1,1)
plot(bins(1:end-1),z)
legend('first half','second half','all trials','Location','northwest')
ylabel('Distance between t and t-1')
title('Jaw Movement')
subplot(2,1,2);
for i=1:length(summary.trial)
    j = ones(1,length(summary.trial(i).centSp))*i;
    h1=scatter(summary.trial(i).centSp,j,6,'filled','r'); %central licks are red
    
     hold on
end

for i=1:length(summary.trial)
    j = ones(1,length(summary.trial(i).LeftSp))*i;
    h2=scatter(summary.trial(i).LeftSp,j,6,'filled','b'); %left licks are blue
     hold on
end
for i=1:length(summary.trial)
    j = ones(1,length(summary.trial(i).RightSp))*i;
    h3=scatter(summary.trial(i).RightSp,j,6,'filled','g'); %right licks are green
     hold on
end
legend([h1,h2,h3],'central','left','right','Location','northwest')
title('Licking')
ylabel('Trial #');
xlabel('Time (sec)')
xlim([-4 6]);
ylim([0 length(summary.trial)]);
suptitle(append(summary.mouseID))

%% licking data
p=1;
tongue_out = find(coords.tongue_p==1)';
%extrace only first frame of a lick for timestamps
temp=flip(tongue_out);
tongue_out_corrected=flip(temp(find(diff((flip(tongue_out)))<-1)));
%kmeans to classify lick direction
t = tongue_out_corrected';

t1(:,1)=coords.tongue_x(t);
t1(:,2)=coords.tongue_y(t);
[idx1, C]=kmeans(t1,3); %kmeans on those coordinates with three centroids, C is centroid and idx is classification of each point in t1
t1(:,3) = idx1;
tongue_out_corrected(2,:) = t1(:,3);
figure;
ax = axes;
hold on
set(ax, 'Ydir', 'reverse')
%t1 contains x and y coordinates for each lick from tongue_out_corrected
%with the third column being the classified cluster from kmeans
e = find(t1(:,3) == 1);
scatter(t1(e,1),t1(e,2)) %plot cluster 1
f = find(t1(:,3) == 2);
scatter(t1(f,1),t1(f,2))%plot cluster 2
g = find(t1(:,3) == 3);
scatter(t1(g,1),t1(g,2)) %plot cluster 3
scatter(C(:,1),C(:,2),'filled') %centroids
title('Kmeans clustering of licking direction')
xlabel('x coordinate')
ylabel('y coordinate')
%labels will change since kmeans chooses clusters randomly
legend('Center','Left','Right','Centroids','Location','Southeast')
cts = [length(e) length(f) length(g)]; %number of licks in each cluster
%extract time stamps of licks from deeplabcut
for i =1:length(summary.trial)
     %idx1 is index of tongue_out_corrected for pulling out catgory from
     %kmeans
     %idx2 is index from the behavior frames to pull out timestamp of licks
     %creating a new field in trial struct that contains lick timestamps
     %and direction as categorized by kmeans clustering
    [~,~,idx1] = intersect( summary.trial(i).Beh_Frame_index,tongue_out_corrected(1,:)); 
    [~,~,idx2] = intersect( tongue_out_corrected(1,:),summary.trial(i).Beh_Frame_index);
    summary.trial(i).lick(1,:) = summary.trial(i).beh_frames(idx2);
    summary.trial(i).lick(2,:) = tongue_out_corrected(2,idx1);
    
end
    
    

%% Plot lick rasters from deep lab cut and detectors
subplot(2,1,1) %plot licks from deeplabcut
for i=1:length(summary.trial)
    p1 = find(summary.trial(i).lick(2,:) == 1); % numbers in p1-p3 correspond to clusters from kmeans
    j1 = ones(1,length(summary.trial(i).lick(1,p1)))*i; %make an array for plotting
    h1=scatter(summary.trial(i).lick(1,p1),j1,6,'filled','b'); %lick raster scatter plot
    p2 = find(summary.trial(i).lick(2,:) == 2);
    j2 = ones(1,length(summary.trial(i).lick(1,p2)))*i;
    h2=scatter(summary.trial(i).lick(1,p2),j2,6,'filled','g');
    p3 = find(summary.trial(i).lick(2,:) == 3);
    j3 = ones(1,length(summary.trial(i).lick(1,p3)))*i;
    h3=scatter(summary.trial(i).lick(1,p3),j3,6,'filled','r');
     hold on
end
title('Licking - DeeplabCut')
ylabel('Trial #');
legend([h1,h2,h3],'central','left','right','Location','Northwest')
xlim([-4 6]);
ylim([0 length(summary.trial)]);
subplot(2,1,2);%plot licks from lick detector
for i=1:length(summary.trial)
    j = ones(1,length(summary.trial(i).centSp))*i;
    h1=scatter(summary.trial(i).centSp,j,6,'filled','r'); %central licks are red
    
     hold on
end

for i=1:length(summary.trial)
    j = ones(1,length(summary.trial(i).LeftSp))*i;
    h2=scatter(summary.trial(i).LeftSp,j,6,'filled','b'); %left licks are blue
     hold on
end
for i=1:length(summary.trial)
    j = ones(1,length(summary.trial(i).RightSp))*i;
    h3=scatter(summary.trial(i).RightSp,j,6,'filled','g'); %right licks are green
     hold on
end
legend([h1,h2,h3],'central','left','right','Location','Northwest')
title('Licking - Detectors')
ylabel('Trial #');
xlim([-4 6]);
ylim([0 length(summary.trial)]);
suptitle(append(summary.mouseID))

%% plot jaw movmemnt along with lick rasters from deeplabcut and detectors
figure;
subplot(3,1,1)
plot(bins(1:end-1),z)
legend('first half','second half','all trials','Location','northwest')
ylabel('Distance between t and t-1')
title('Jaw Movement')
subplot(3,1,2) %plot licks from deeplabcut
for i=1:length(summary.trial)
    p1 = find(summary.trial(i).lick(2,:) == 1); % numbers in p1-p3 correspond to clusters from kmeans
    j1 = ones(1,length(summary.trial(i).lick(1,p1)))*i; %make an array for plotting
    h1=scatter(summary.trial(i).lick(1,p1),j1,6,'filled','r'); %lick raster scatter plot
    p2 = find(summary.trial(i).lick(2,:) == 2);
    j2 = ones(1,length(summary.trial(i).lick(1,p2)))*i;
    h2=scatter(summary.trial(i).lick(1,p2),j2,6,'filled','b');
    p3 = find(summary.trial(i).lick(2,:) == 3);
    j3 = ones(1,length(summary.trial(i).lick(1,p3)))*i;
    h3=scatter(summary.trial(i).lick(1,p3),j3,6,'filled','g');
     hold on
end
title('Licking - DeeplabCut')
ylabel('Trial #');
legend([h1,h2,h3],'central','left','right','Location','Northwest')
xlim([-4 6]);
ylim([0 length(summary.trial)]);
subplot(3,1,3);%plot licks from lick detector
for i=1:length(summary.trial)
    j = ones(1,length(summary.trial(i).centSp))*i;
    h1=scatter(summary.trial(i).centSp,j,6,'filled','r'); %central licks are red
    
     hold on
end

for i=1:length(summary.trial)
    j = ones(1,length(summary.trial(i).LeftSp))*i;
    h2=scatter(summary.trial(i).LeftSp,j,6,'filled','b'); %left licks are blue
     hold on
end
for i=1:length(summary.trial)
    j = ones(1,length(summary.trial(i).RightSp))*i;
    h3=scatter(summary.trial(i).RightSp,j,6,'filled','g'); %right licks are green
     hold on
end
legend([h1,h2,h3],'central','left','right','Location','northwest')
title('Licking - Detectors')
ylabel('Trial #');
xlim([-4 6]);
ylim([0 length(summary.trial)]);
suptitle(append(summary.mouseID))


