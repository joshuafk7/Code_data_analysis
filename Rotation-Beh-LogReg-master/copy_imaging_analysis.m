clear all;clc

% Imaging analysis for GC project
%% rigister your image
 %edit regtiff
%% load imaging data and ROI
%load('D:\Behavior\Discrimination\Imaging\JK105\200126\neuron.mat')
% neuron = results;
load('/Users/ayeshav/Documents/Rotation3/02-Apr_16_30_34.mat')

%neuron = neuron_new2;
clearvars -except neuron
%% load the event
file = dir('*.rhd');
a = {'Suc','S_75_25','S_55_45','S_45_55','S_25_75', 'NaCl'};
% a = {'Suc','S_85_15','S_65_35','S_35_65','S_15_85', 'NaCl'};
b = [2 2 2 1 1 1];
[data,trial,summary] = process_intan_v4_behavior_only_edit(file.name, a, b);


%% get inferredSpikeTimes

spikeTimes = getSpikeTimes(data.imaging_frames, full(neuron.S));


%% spatial mapping
%[correlationmatrix,center,sortedD] = Spatial_organization_CNMFE_Josh_copy(neuron);
%h= heatmap(correlationmatrix);
%h.Colormap = parula;

%% binned data - create 4d matrix 
%4D array called binnedC, dim1=trial, dim2=bin, dim3=taste, dim4=neuron

win_sample = [-2 3]; binWidth = 0.121;
win_action = [-3 2]; 

spikeData = [];

for i = 1:size(neuron.S,1)

spikeData_sampleTime(:,:,i) = getinferredSpikeAct(spikeTimes(i), [trial.sampleTime], win_sample, binWidth, []);
spikeData_actionTime(:,:,i) = getinferredSpikeAct(spikeTimes(i), [trial.actionOnset], win_action, binWidth, []);

end


x = -2:binWidth:3-binWidth;
 
var1 = nanmean(nanmean(spikeData_sampleTime(([trial.tasteStim]==1),:,:)./binWidth,1),3);
var2 = nanmean(nanmean(spikeData_sampleTime(([trial.tasteStim]==7),:,:)./binWidth,1),3);


plot(x,smoothdata(var1./binWidth,'gaussian',10),'r-','LineWidth',2);hold on;
plot(x,smoothdata(var2./binWidth,'gaussian',10),'b-','LineWidth',2);hold on;


plot([0 0],ylim(),'k--')

xlabel('time from sample onset (s)')
