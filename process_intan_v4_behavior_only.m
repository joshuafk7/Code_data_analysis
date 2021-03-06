%% process data recorded from Intan board
%input variables
    %filename - the .rhd file output from the intan, one recorded for each
    %behavioral session
    %excel_tastes - tastes in each line example {'Suc', 'S_75_25','S_55_45','S_45_55','S_25_75','NaCl'}
    %excel_directions - correct direction for each line, 1 left, 2 right,
    % example [1 1 1 2 2 2]
%output
    %all other outputs are saved in the summary variable but you can also output
    %things individually if you want
function [tastes,unit, data,trial, summary] = process_intan_v4_behavior_only(file,excel_tastes,excel_directions)
% file = dir('RVKC441_191001_112219.rhd');
dataRaw = read_Intan(file.name);
% data = read_intan_batch;
%% extract the event data
thr = .5;
[data.centSp(1,:),data.centSp(2,:)] = Timing_onset_offset(dataRaw.analog(1,:), dataRaw.ts, thr,30,0); % get the central licks
[data.LeftSp(1,:),data.LeftSp(2,:)] = Timing_onset_offset(dataRaw.analog(2,:), dataRaw.ts, thr,30,0);
[data.RightSp(1,:),data.RightSp(2,:)] = Timing_onset_offset(dataRaw.analog(3,:), dataRaw.ts, thr,30,0);
[data.beh_frames(1,:),data.beh_frames(2,:)] = Timing_onset_offset(dataRaw.analog(4,:), dataRaw.ts, thr,20,0);

[data.Forward(1,:),data.Forward(2,:)]              = Timing_onset_offset(dataRaw.event(10,:), dataRaw.ts, 0.5,3000,0);
[data.Backward(1,:),data.Backward(2,:)]              = Timing_onset_offset(dataRaw.event(11,:), dataRaw.ts, 0.5,30,0);
[data.Up(1,:),data.Up(2,:)]                        = Timing_onset_offset(dataRaw.event(12,:), dataRaw.ts, 0.5,30,0);
[data.Down(1,:),data.Down(2,:)]                    = Timing_onset_offset(dataRaw.event(13,:), dataRaw.ts, 0.5,30,0);

%%
%% events
A = cd;
B=file.name;
summary.mouseID  = B(1:5);
summary.date     = A(end-5:end);
clear A;
% %%
% [data.T_1(1,:),data.T_1_off(2,:)]     = Timing_onset_offset(dataRaw.event(1,:), dataRaw.ts, 0.5,30,0);
% [data.T_2,~]     = Timing_onset_offset(dataRaw.event(2,:), dataRaw.ts, 0.5,30,0);
% [data.T_3,~]     = Timing_onset_offset(dataRaw.event(3,:), dataRaw.ts, 0.5,30,0);
% [data.T_4,~]     = Timing_onset_offset(dataRaw.event(4,:), dataRaw.ts, 0.5,30,0);
% [data.T_5,~]     = Timing_onset_offset(dataRaw.event(5,:), dataRaw.ts, 0.5,30,0);
% [data.T_6,~]     = Timing_onset_offset(dataRaw.event(6,:), dataRaw.ts, 0.5,30,0);
% [data.T_7,~]     = Timing_onset_offset(dataRaw.event(7,:), dataRaw.ts, 0.5,30,0);
for i = 1:7
    if ~isempty(Timing_onset_offset(dataRaw.event(i,:), dataRaw.ts, 0.5,30,0))
    [data.(append('T_',num2str(i)))(1,:),data.(append('T_',num2str(i)))(2,:)] = ...
        Timing_onset_offset(dataRaw.event(i,:), dataRaw.ts, 0.5,30,0);
    end
end
data.T_4 = [];
data = rmfield(data, 'T_4'); %remove rinses
[data.R_1(1,:), data.R_1(2,:)]            = Timing_onset_offset(dataRaw.event(8,:), dataRaw.ts, 0.5,30,0);

[data.L_1(1,:), data.L_1(2,:)]            = Timing_onset_offset(dataRaw.event(9,:), dataRaw.ts, 0.5,30,0);   
if ~isempty(find(dataRaw.event(14,:)))
[data.imaging_frames_raw(1,:), data.imaging_frames_raw(2,:)]            = Timing_onset_offset(dataRaw.event(14,:), dataRaw.ts, 0.5,30,0);   
data.imaging_frames = data.imaging_frames_raw(:,1:4:end);
end


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
     if isempty(data.(data_signals{i}))  %remove empty fields
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
 assigned_tastes = excel_tastes;
 assigned_direction = excel_directions;
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
start=-4000; %how long before event to start trial
finish =6000; %how long after event to start trial


for i =1:length(data_signals)
    if ~isempty(data.(data_signals{i}))
    unit.(data_signals{i}) = spike2eventRasteandPSTH_NP_Josh(data.(data_signals{i}),data.firstLick,100,start,finish);
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
%             if idx ==1
%                trial(idx).ITI = 0;
%             else
%                 trial(idx).ITI = data.firstLick(i)-data.Down(i-1);
%             end
        end
        for a = 1:length(tastes)
            if ~isempty(unit.(tastes{a}).spikeraster(i).times)
                trial(idx).TasteID = tastes{a}; %extract taste ID for each trial
            end
        end
        for b = 1:length(tastes)
            if isfield(unit, 'beh_frames')==1
            trial(idx).beh_frames_index = unit.beh_frames.spikeraster(i).index; %extrac frame index for alignment with imaging
            end
            if isfield(unit,'imaging_frames')==1
            trial(idx).imaging_frames_index = unit.imaging_frames.spikeraster(i).index; %extrac frame index for alignment with imaging
            end
        end
        idx=idx+1;
%     end
    
end
%% list trials with no lateral licks 

j=1;
lateralmiss=[];
for i=1:length(trial)
    
    A = isempty(trial(i).LeftSp);
    B = isempty(trial(i).RightSp);
    
    if A&&B
        lateralmiss(j) = i;
        j=j+1;
    end
end

%% problem trials (random event recorded by national instruments results in event in both left and right water delivery
j=1;
problemTrial = [];
for i =1:length(trial)
    if ~isempty(trial(i).R_1) & ~isempty(trial(i).L_1)
        problemTrial(j) = i;
        j=j+1;
        
    end
end

%% remove trial where there is no taste ID

j=1;
noTasteID = []; 
for i =1:length(trial)
    if isempty(trial(i).TasteID)
        noTasteID(j) = i;
        j=j+1; 
    
    end
end
%% remove these trials from trial struct

a=[problemTrial  lateralmiss noTasteID];

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
    trial(i).bias = Lerror/Lcount - Rerror/Rcount;
end
summary.bias = Lerror/Lcount - Rerror/Rcount;
    

%% summary performance for each taste
T_correct = zeros(1,length(tastes));
T_count = zeros(1,length(tastes));
t_performance = zeros(1,length(tastes));
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
% for j=1:length(tastes)
%     summary.(append((tastes{j}),'_performance')) = T_correct(j)/T_count(j);
% end
t_performance = T_correct./T_count;
summary.ind_performance = t_performance;


%% summary total performance
%all outputs of this function are also saved in the summary struct
summary.total_performance = trial(length(trial)).total_performance;
summary.numTrials = length(trial);
summary.tastes = tastes;
summary.directions = excel_directions;
summary.trial = trial;
summary.unit = unit;
summary.data = data;
summary.problemTrial = problemTrial;
summary.lateralmiss = lateralmiss;
summary.noTasteID = noTasteID;
    
end