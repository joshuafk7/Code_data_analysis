function [unit] = spike2eventRasteandPSTH_NP (Spike, Event, Binsize, Pre, Post)
%function [unit] = spike2eventhistraster (Spike, Event, Binsize, Pre, Post)
% Haixin Liu 2014
% this function is used to extract info for raster and psth for a given
% spike and event in a single session
% INPUT:
%      Spike: vector, timestamps of all spikes of a unit
%      Event: vector, timestamps of all spikes of a event
%      Binsize: in ms, 0 point for raster and psth
%      Pre, Post: in ms for time peirod prior and after Event
% OUTPUT: a structure
%      unit.
%           Spike: record input Spike
%           Event: record input Event
%           Binsize: record input Binsize
%           timebin: record bins used to bin raster plot
%           sc_sum: spike count per bin, sum along trial
%           scmatrix: spike count matrix
%           FRmatrix: firing rate matrix
%           sc_avg: spike count average
%           FR_avg: firing rate average
%           timepoint: timepoint for every bin (center)
%           spikeraster: used for raster plots; spikeraster(trial#).time :
%                        a vector of timestamps for this trial aligned to
%                        event timestamp
%           
%%
unit = [];
unit.Spike = Spike;
unit.Event = Event;
unit.Binsize = Binsize;

pre = Pre*0.001;
post = Post*0.001;
bin = Binsize*0.001;
unit.timebin = (pre:bin:post);
edge = unit.timebin';
%unit.timebin can serve as edge for histc
n = zeros(length(unit.Event), length(edge));
% n is a matrix to record histc result of every trial
for i = 1:length(unit.Event)
    j = 1; k = 0;
    spike =[];
    %spike as storing vector to store every trial spike timestamp to avoid
    %0 which generated when putting trials in matrix unit.spike; reset
    %every trial
    %reset j
    for j=1:length(unit.Spike)
        if (unit.Spike(j) >= unit.Event(i)+pre)&&(unit.Spike(j) <= unit.Event(i)+post)
            k=k+1;
            spike(k)= unit.Spike (j)- unit.Event(i);
            %unit.spike (i, k)= unit.Spike (j)- unit.Event(i);
        else
        end
    end
    if (isempty(spike))
        n(i,:)=0;
    else
        n(i,:) = histc (spike, edge);
    end
    %count numbers falling in every bin
    %The last bin counts any values of x that match edges(end)
end
n(:,end) = [];
unit.sc_sum = sum(n);
unit.scmatrix = n;
unit.FRmatrix = unit.scmatrix./bin;
unit.sc_avg = unit.sc_sum./length(unit.Event);
unit.FR_avg = unit.sc_avg./bin;
unit.timepoint = unit.timebin+bin/2;%because every bin starts and ends, so make the middle point as timepoint of the bin
unit.timepoint(:,end) = [];



spikeraster = [];

%raster part
for i = 1:length(unit.Event)
    j = 1;
    spikeraster(i).times = [];
    %reset j
    k=1;
    f=1;
    
    for j=1:length(unit.Spike)
        if (unit.Spike(j) >= unit.Event(i)+pre)&&(unit.Spike(j) <= unit.Event(i)+post)
            
            spikeraster(i).times(k)= unit.Spike (j)- unit.Event(i); %time stamp of each spike relative to event
            spikeraster(i).index(f) = j; %josh added to get index of frames for imaging
            f=f+1; %keep track of indicies
            k=k+1; %keep track of times
            %unit.spike (i, k)= unit.Spike (j)- unit.Event(i);
        else
        end
    end
end
unit.spikeraster = spikeraster;
% if length(unit.FR_avg) == 1
%     unit.preFR = 0;
%     unit.overallFR = 0;
% else
% unit.preFR = mean(unit.FR_avg(unit.timepoint<=0));
% %fprintf (['prestimulus FR is ',num2str(unit.preFR),' Hz \n'])
% unit.overallFR = mean(unit.FR_avg);
% end
end
%fprintf (['overall FR is ',num2str(unit.overallFR),' Hz \n'])

