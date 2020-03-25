%% process data recorded from Intan board
function [tastes,unit, data,trial, summary] = process_intan_v2_josh(filename);
% filename = 'D:\Behavior\Discrimination\JK062\190914\JK062_190914_161300.rhd';
file = dir('*.rhd');
dataRaw = read_Intan(file.name);
% data = read_intan_batch;
%% extract the event data
thr = .5;
[data.centSp(1,:),data.centSp(2,:)] = Timing_onset_offset(dataRaw.analog(1,:), dataRaw.ts, thr,30,0); % get the central licks
[data.LeftSp(1,:),data.LeftSp(2,:)] = Timing_onset_offset(dataRaw.analog(2,:), dataRaw.ts, thr,30,0);
[data.RightSp(1,:),data.RightSp(2,:)] = Timing_onset_offset(dataRaw.analog(3,:), dataRaw.ts, thr,30,0);
% [data.Forward,~]              = Timing_onset_offset(dataRaw.event(10,:), dataRaw.ts, 0.5,3000,0);
% [data.Backward,~]              = Timing_onset_offset(dataRaw.event(11,:), dataRaw.ts, 0.5,30,0);
% [data.Up,~]                        = Timing_onset_offset(dataRaw.event(12,:), dataRaw.ts, 0.5,30,0);
% [data.Down,~]                    = Timing_onset_offset(dataRaw.event(13,:), dataRaw.ts, 0.5,30,0);
[data.Frames(1,:),data.Frames(2,:)]                    = Timing_onset_offset(dataRaw.event(14,:), dataRaw.ts, 0.5,30,0);

%% mouse and date info

% %% events
 A = cd;
summary.mouseID  = A(28:32);
summary.date     = A(34:end);
% clear A;
%%
for i = 1:7
    if ~isempty(Timing_onset_offset(dataRaw.event(i,:), dataRaw.ts, 0.5,30,0))
    [data.(append('T_',num2str(i)))(1,:),data.(append('T_',num2str(i)))(2,:)] = ...
        Timing_onset_offset(dataRaw.event(i,:), dataRaw.ts, 0.5,30,0);
    end
end
[data.R_1(1,:), data.R_1(2,:)]            = Timing_onset_offset(dataRaw.event(8,:), dataRaw.ts, 0.5,3000,0);
[data.L_1(1,:), data.L_1(2,:)]            = Timing_onset_offset(dataRaw.event(9,:), dataRaw.ts, 0.5,3000,0);        
% [data.T_1,data.T_1_off]     = Timing_onset_offset(dataRaw.event(1,:), dataRaw.ts, 0.5,30,0);
% [data.T_2,~]     = Timing_onset_offset(dataRaw.event(2,:), dataRaw.ts, 0.5,30,0);
% [data.T_3,~]     = Timing_onset_offset(dataRaw.event(3,:), dataRaw.ts, 0.5,30,0);
% [data.T_4,~]     = Timing_onset_offset(dataRaw.event(4,:), dataRaw.ts, 0.5,30,0);
% [data.T_5,~]     = Timing_onset_offset(dataRaw.event(5,:), dataRaw.ts, 0.5,30,0);
% [data.T_6,~]     = Timing_onset_offset(dataRaw.event(6,:), dataRaw.ts, 0.5,30,0);
% [data.T_7,~]     = Timing_onset_offset(dataRaw.event(7,:), dataRaw.ts, 0.5,30,0);


        
% data_signals = {'centSp','LeftSp','RightSp','Forward','Backward','Up','Down','Frames','T_1','T_2','T_3','T_4','T_5','T_6','T_7','R_1','L_1'};
%% remove NI errors
names = fieldnames(data);
for i = 1:length(names)
    data.(names{i})(3,:)= data.(names{i})(1,:)-data.(names{i})(2,:);
end
for i =1:length(names)
   A=find(data.(names{i})(3,:) == 0);
   data.(names{i})(:,A) = [];
end
for i = 1:length(names)
   data.(names{i})(2:3,:) = []; 
end
%% find trials with no central sampling 
%only need these if you are aligning to something other than lick signal
% q = spike2eventRasteandPSTH_NP_Josh(data.centSp,data.Forward,100,-2000,2000); 
% j=1;
% Omit_trials=[];
% for i =1:length(q.spikeraster)
%     if  isempty(q.spikeraster(i).times) | q.spikeraster(i).times ==0 | q.spikeraster(i).times(1)<0 %if there was a forward with no lick afterwards
%         Omit_trials(j) = i;
%         j=j+1;
%     end
% end
    
%% Frames adjusted takes into account downsampling of video
data.Frames_adjusted = data.Frames(1,1:end);
data.Frames = [];
%% this section extracts first lick of each trial for alignment purposes
%it also means that trials without central sp licks are not extracted
x=2;
data.firstLick(1) = data.centSp(1);
for i =2:length(data.centSp) 
    if data.centSp(i)-data.centSp(i-1) <0.5
        continue
    else
        data.firstLick(x) = data.centSp(i);
        x=x+1;
    end
end

%% for

%% remove empty tastes
 data_signals = fieldnames(data);
 b=1;
 for i=1:length(data_signals)
     if isempty(data.(data_signals{i})) && ~isequal(data_signals{i},'Frames') %remove empty fields
        data= rmfield(data, (data_signals{i}));
     end
 end
 data_signals = fieldnames(data);
 
 %% extract tastes from data struc and allow renaming
 for i =1:length(data_signals)
     if  startsWith(data_signals(i) ,'T')
         tastes(b) = data_signals(i);
         b=b+1;
     end
 end
 tastes = sort(tastes);
 assigned_tastes = input('Enter tastants in order from line 1 \n');
 assigned_direction = input('Input direction of tastes in array (1 is left, 2 is right)\n');
 for i =1:length(assigned_tastes)
     data.(assigned_tastes{i}) = data.(tastes{i});
 end
 
 %% remove tastes with T label, now tastes are only labeled with user assigned names
 for i =1:length(tastes)
     data= rmfield(data, (tastes{i}));
 end
data_signals = fieldnames(data);   
tastes = assigned_tastes;
clear -except assigned_tastes;
 
%% unit struct contains spike rasters for all signals aligned to first lick of central spout
start=-2000; %how long before event to start trial
finish =4500; %how long after event to start trial

startFrames = -5000; %get frames from longer period of time, this is is to prevent issues with trial misalignment
finishFrames = 10000;

for i =1:length(data_signals)
    if ~isequal(data_signals{i},'Frames_adjusted')
        unit.(data_signals{i}) = spike2eventRasteandPSTH_NP_Josh(data.(data_signals{i}),data.firstLick,100,start,finish);
    else
%     if ~isempty(data.Frames)
        unit.Frames_adjusted = spike2eventRasteandPSTH_NP_Josh(data.Frames_adjusted,data.firstLick,100,startFrames,finishFrames); %grab froms from longer period of time
%     end
    end
end
%% Reorganize in trial structure
idx=1;
for i = 1:length(data.firstLick)
    %if statement gets rid of trials where not all frames were recorded by
    %15 refers to approx framerate
    %start and finish from above divided by 1000 to get trial length in sec
    %frame rate shouldn't vary by more than 10 less than trial length*hz
    %unless frames were dropped
%     if length(unit.Frames.spikeraster(i).times) > (((-start +finish) /1000)*15)-10 
        trial(idx).firstLick = data.firstLick(i);
        for j = 1:length(data_signals)
            trial(idx).(data_signals{j}) = unit.(data_signals{j}).spikeraster(i).times;
            trial(idx).OverallFirstLick = data.firstLick(i);
        end
        for b = 1:length(data_signals)
            trial(idx).Frame_index = unit.Frames_adjusted.spikeraster(i).index; %extrac frame index for alignment with imaging
        end
        for a = 1:length(tastes)
            if ~isempty(unit.(tastes{a}).spikeraster(i).times)
                trial(idx).TasteID = tastes{a}; %extract taste ID for each trial
            end
        end
        idx=idx+1;
%     end
    
end

%delete first and last trials (they dont have the right number of frames
trial(1) = [];
trial(length(trial)) = [];

%% make number of frames and frame index uniform for all trials (this required deleting one frame from approx half trials

% find(min(length(trial.Frames_adjusted)))
%% list trials with no lateral licks 

j=1;
summary.lateralmiss=[];
for i=1:length(trial)
    
    A = isempty(trial(i).LeftSp);
    B = isempty(trial(i).RightSp);
    
    if A&&B
        summary.lateralmiss(j) = i;
        j=j+1;
    end
end

%% problem trials (random event recorded by national instruments results in event in both left and right water delivery
j=1;
summary.problemTrial = [];
for i =1:length(trial)
    if ~isempty(trial(i).R_1) & ~isempty(trial(i).L_1)
        summary.problemTrial(j) = i;
        j=j+1; 
    end
end

%% remove trial where there is no taste ID
summary.noTasteID = []; 
j=1;
for i =1:length(trial)
    if isempty(trial(i).TasteID)
        summary.noTasteID(j) = i;
        j=j+1; 
    end
end
%% remove these trials from trial struct
a=[summary.problemTrial  summary.lateralmiss summary.noTasteID];

trial(a) = [];

%% add if it was a left or right trial based on user assigned directions

for i=1:length(trial)
    for j=1:length(tastes)
        if convertCharsToStrings(trial(i).TasteID(:)) == tastes{j}
            trial(i).L_R_trial = assigned_direction(j);
        end
    end
end

%% add correct or incorrect trials based on if water was delivered at lateral spout

for i =1:length(trial)
    if isempty(trial(i).L_1) & isempty(trial(i).R_1)
        trial(i).correct_choice = 0;
    else
        trial(i).correct_choice = 1;
    end
end


%% total performance
x=0;
for i=1:length(trial)
    if trial(i).correct_choice == 1
        x=x+1;
    end
    trial(i).total_performance = x/i;
end
%% performance for each taste
x=zeros(1,length(tastes));
trialcount=zeros(1,length(tastes));
b=1;
for j=1:length(tastes)
    for i=1:length(trial)
        if trial(i).correct_choice == 1 & convertCharsToStrings(trial(i).TasteID(:)) == tastes{j}
            x(b)=x(b)+1;   
        end
        
        if convertCharsToStrings(trial(i).TasteID(:)) == tastes{j}
            trialcount(b)=trialcount(b)+1;
            trial(i).((append((tastes{j}),'_performance'))) = x(b)/trialcount(b);
            
        else
            trial(i).((append((tastes{j}),'_performance'))) = [];
        end
    end
    b=b+1;
end

%% bias calculation are they making more errors in one direction than the other

Lerror = 0;
Rerror=0;
Lcount = 0;
Rcount = 0;
for i =1:length(trial)
    if trial(i).correct_choice == 0 & trial(i).L_R_trial == 1
        Lerror = Lerror+1;
    end
    if trial(i).correct_choice == 0 & trial(i).L_R_trial == 2
        Rerror = Rerror+1;
    end
    if trial(i).L_R_trial == 1
        Lcount = Lcount+1;
    end
    if trial(i).L_R_trial == 2
        Rcount = Rcount+1;
    end
end
summary.bias = Lerror/Lcount - Rerror/Rcount;
    

%% summary performance for each taste
T_correct = zeros(1,length(tastes));
T_count = zeros(1,length(tastes));
for j=1:length(tastes)
    for i =1:length(trial)
        if convertCharsToStrings(trial(i).TasteID) == tastes{j}
            T_count(j) = T_count(j)+1;
            if trial(i).correct_choice == 1
                T_correct(j) = T_correct(j)+1;
            end
        end   
    end
end
for j=1:length(tastes)
    summary.(append((tastes{j}),'_performance')) = T_correct(j)/T_count(j);
end



%% summary total performance
summary.total_performance = trial(length(trial)).total_performance;
    
end