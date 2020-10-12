% figure(9987)
% hold on
ROItest = ReadImageJROI('RoiSet_interneurons.zip');
% load('neuron_filtered.mat');
% neuron=neuron_new;
[dim1, dim2 ] = size(imread('JK146_mcherry_sameplane.jpeg'));
for i=1:length(ROItest)
   A = length(ROItest{1,i}.mnCoordinates );
   ROItest{1,i}.mnCoordinates(A+1,:)=ROItest{1,i}.mnCoordinates(1,:);
end
i=1;
for i=1:length(ROItest)
mcherry_ROIs(:,:,i) = roipoly(dim1,dim2,ROItest{1,i}.mnCoordinates(:,1),ROItest{1,i}.mnCoordinates(:,2));

end
Coor = neuron.show_contours(0.8);
for i=1:length(neuron.Coor)
CNMF_ROI1(:,:,i) = roipoly(dim1,dim2,neuron.Coor{i,1}(1,:),neuron.Coor{i,1}(2,:));

end
%%
A = neuron.reshape(neuron.A, 2);
B = full(neuron.A);
CNMF_ROIs = reshape(B , size(neuron.Cn,1), size(neuron.Cn,2), size(B,2));



ROI1_compressed=NaN(1000,size(mcherry_ROIs,3));
for i=1:size(mcherry_ROIs,3)
temp1 = find(mcherry_ROIs(:,:,i));
temp2=1000-length(temp1);
ROI1_compressed(:,i)=[temp1; NaN(temp2,1)];
end

% length(find(CNMF_ROIs(:,:,1)))
% length(find(CNMF_ROIs(:,:,1)>2))


% 
% ROI2_compressed=NaN(1000,size(CNMF_ROIs,3));
% for i=1:size(CNMF_ROIs,3)
% temp1 = find(CNMF_ROIs(:,:,i)>2);
% temp2=1000-length(temp1);
% ROI2_compressed(:,i)=[temp1; NaN(temp2,1)];
% end


ROI3_compressed=NaN(1000,size(CNMF_ROI1,3));
for i=1:size(CNMF_ROI1,3)
temp1 = find(CNMF_ROI1(:,:,i));
temp2=1000-length(temp1);
ROI3_compressed(:,i)=[temp1; NaN(temp2,1)];
end

p=1;
for i=1:size(mcherry_ROIs,3)
    for j=1:size(CNMF_ROI1,3)
        if intersect(ROI1_compressed(:,i),ROI3_compressed(:,j))~=0
            overlapping_pair(p,1)=i;
            overlapping_pair(p,2)=j;
            
            p=p+1;
        end
    end
end
for i=1:length(overlapping_pair)
    overlapping_pair(i,3)=sum(~isnan(ROI1_compressed(:,overlapping_pair(i,1))));
    overlapping_pair(i,4)=sum(~isnan(ROI3_compressed(:,overlapping_pair(i,2))));
    overlapping_pair(i,5) = length(intersect(ROI1_compressed(:,overlapping_pair(i,1)),ROI3_compressed(:,overlapping_pair(i,2))));
end
overlapping_pair(:,6)=overlapping_pair(:,5)./overlapping_pair(:,3);

GC=[];
[GC(:,2),GC(:,1)]=groupcounts(overlapping_pair(:,1));

unique_interneuron_idx = find(GC(:,2)==1);
unique_interneuron = GC(unique_interneuron_idx,1);
all_neurons = 1:length(neurons);
excitatory = setdiff(all_neurons, overlapping_pair(:,2));
%% 
p=1;
for i=1:length(overlapping_pair)
if ismember(overlapping_pair(i,1),unique_interneuron)
    overlapping_pair_filtered(p,:) = overlapping_pair(i,:);
    p=p+1;

end
end


p=1;
overlapping_pair_filtered2=[];
for i=1:length(overlapping_pair_filtered)
if overlapping_pair_filtered(i,6)>=.75 % set overlap threshold
    overlapping_pair_filtered2(p,:) = overlapping_pair_filtered(i,:);
    p=p+1;

end
end

%% 
overlapping_pair_filtered3=[];
sample_thresholds=[]; 
num_neurons_threshold=[];
sample_thresholds = linspace(0,1,21);
for j=1:(length(sample_thresholds))-5
    p=1;
for i=1:length(overlapping_pair_filtered)
if overlapping_pair_filtered(i,6)>=sample_thresholds(j) && overlapping_pair_filtered(i,6)<=sample_thresholds(j+5) % set overlap threshold
    overlapping_pair_filtered3{1,j}(p,:) = overlapping_pair_filtered(i,:);
    p=p+1;

end
end
num_neurons_threshold(j) = length(overlapping_pair_filtered3{1,j});
end

figure(837)
bar(sample_thresholds(1:end-5),num_neurons_threshold)
xlabel('Amount Overlap')
ylabel('Num Interneurons Identified')
title('Testing threshold for Interneuron Inclusion')
 %%
 interneurons=[];
 task_related_interneurons=[];
%  overlapping_pair_filtered3 = overlapping_pair_filtered2;
 for i=1:length(overlapping_pair_filtered3)
 interneurons{1,i} = overlapping_pair_filtered3{1,i}(:,2);
 end
 
 for i=1:length(overlapping_pair_filtered3)
     temp111=intersect(interneurons{1,i},taste_resp_neurons);
    temp222=intersect(interneurons{1,i},delay_resp_neurons);
    temp333=intersect(interneurons{1,i},choice_resp_neurons);
    task_related_interneurons{1,i} = unique([temp111;temp222;temp333]);
 end
 %%
 figure(790)
for i=1:length(task_related_interneurons)
subplot(4,4,i)
popavgtaste1 = mean(squeeze(nanmean(plotbinnedC7(:,:,2,task_related_interneurons{1,i}),1)),2)';
popavgtaste2 = mean(squeeze(nanmean(plotbinnedC7(:,:,2,task_related_neurons_all),1)),2)';
plot(plotbins(1:end-1),popavgtaste1);
ylim = [0 .1];

hold on
plot(plotbins(1:end-1),popavgtaste2);
y=ylim;
x=xlim;
rectangle('Position',[0 0 2 y(2)],'EdgeColor', 'none', 'FaceColor', [.8, .8, .8, 0.5])
line([0 0],[0 y(2)],'Color','r','LineStyle','--')
title(append('Threshold =',num2str(sample_thresholds(i)),' - ',num2str(sample_thresholds(i+5))))
text(x(1)-4, y(2)-.02,num2str(length(task_related_interneurons{1,i})))
hold off
% 

% if i==1
% legend('Interneurons','All Neurons','Location','northwest')
% end
% 
if mod(i,4) == 1
    ylabel('Mean Spiking')
end
if i>=13
    xlabel('Time (sec)')
end
end
sgtitle('Threshold Ranges - Interneuron Identification')
%% averaging for statistics
diff_spiking=[];
pval_thresholds=[];
for i = 1:length(task_related_interneurons)
trial_avg_interneuron = nanmean(squeeze(binnedC7(:,3,2,task_related_interneurons{1,i})),2);
trial_avg_all_neurons = nanmean(squeeze(binnedC7(:,3,2,task_related_neurons_all)),2);
diff_spiking(i) = nanmean(trial_avg_interneuron)-nanmean(trial_avg_all_neurons);
pval_thresholds(i)=ranksum(trial_avg_interneuron,trial_avg_all_neurons);
end

%%
figure(555)
for i=1:length(task_related_interneurons)
    
    
    for j=1:1000
        
        temp1 = randperm(length(task_related_neurons_all),num_neurons_threshold(i)) ;
        trial_avg_interneuron = nanmean(squeeze(binnedC7(:,3,2,task_related_interneurons{1,i})),2);
        trial_avg_all_neurons = nanmean(squeeze(binnedC7(:,3,2,task_related_neurons_all(temp1))),2);
        diff_spiking_histo(j) = nanmean(trial_avg_interneuron)-nanmean(trial_avg_all_neurons);
        
    end
    stdev2x_plus_mean = mean(diff_spiking_histo)+std(diff_spiking_histo)*2;
    stdev2x_minus_mean = mean(diff_spiking_histo)-std(diff_spiking_histo)*2;
    subplot(4,4,i)
    histogram(diff_spiking_histo,'Normalization','probability')
    ylim =[0 .2];
    xlim = [ -.2  .2];
    y=ylim;
    x=xlim;
    line([stdev2x_plus_mean stdev2x_plus_mean],[0 y(2)],'Color','r','LineStyle','--')
    line([stdev2x_minus_mean stdev2x_minus_mean],[0 y(2)],'Color','r','LineStyle','--')
    line([diff_spiking(i) diff_spiking(i)],[0 y(2)],'Color','g','LineStyle','--')
%     title(append('Threshold = ',num2str(sample_thresholds(i)),' - ',num2str(sample_thresholds(i+1))))
    title(append(num2str(sample_thresholds(i)),' - ',num2str(sample_thresholds(i+5))))
%     text(x(1), y(2),num2str(length(task_related_interneurons{1,i})))
    length(task_related_interneurons{1,i});
end
sgtitle('Threshold Ranges')
%%
task_related_neurons_all = unique([taste_resp_neurons;delay_resp_neurons;choice_resp_neurons]);




%%
CNMF_compressed =zeros(size(CNMF_ROIs,1),size(CNMF_ROIs,2));
for i=1:size(overlapping_pair_filtered2,1)
   CNMF_compressed(:,:) = CNMF_compressed(:,:) +CNMF_ROIs(:,:,overlapping_pair_filtered2(i,2));
end

figure(1547)
imagesc(CNMF_compressed)

mcherry_compressed =zeros(size(mcherry_ROIs,1),size(mcherry_ROIs,2));
for i=1:size(overlapping_pair_filtered2,1)
   mcherry_compressed(:,:) = mcherry_compressed(:,:) +mcherry_ROIs(:,:,overlapping_pair_filtered2(i,1));
end

figure(8236)
imagesc(mcherry_compressed)
