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
function [data,trial, summary] = process_intan_v4_behavior_only_edit(filename,excel_tastes,excel_directions)
% filename = 'D:\Behavior\Discrimination\JK062\190914\JK062_190914_161300.rhd';
file = dir('*.rhd');
dataRaw = read_Intan(file.name);
% data = read_intan_batch;
%% extract the event data
thr = .5;
[data.centSp(1,:),data.centSp(2,:)] = Timing_onset_offset(dataRaw.analog(1,:), dataRaw.ts, thr,30,0); % get the central licks
[data.LeftSp(1,:),data.LeftSp(2,:)] = Timing_onset_offset(dataRaw.analog(2,:), dataRaw.ts, thr,30,0);
[data.RightSp(1,:),data.RightSp(2,:)] = Timing_onset_offset(dataRaw.analog(3,:), dataRaw.ts, thr,30,0);
[data.beh_frames(1,:),data.beh_frames(2,:)] = Timing_onset_offset(dataRaw.analog(4,:), dataRaw.ts, thr,20,0);


%%
%% events
A = cd;
B=file.name;
summary.mouseID  = B(1:5);
summary.date     = A(end-5:end);
clear A;

for i = 1:7
    if ~isempty(Timing_onset_offset(dataRaw.event(i,:), dataRaw.ts, 0.5,30,0))
    [data.(append('T_',num2str(i)))(1,:),data.(append('T_',num2str(i)))(2,:)] = ...
        Timing_onset_offset(dataRaw.event(i,:), dataRaw.ts, 0.5,30,0);
    end
end

[data.R_1(1,:), data.R_1(2,:)] = Timing_onset_offset(dataRaw.event(8,:), dataRaw.ts, 0.5,30,0);

[data.L_1(1,:), data.L_1(2,:)] = Timing_onset_offset(dataRaw.event(9,:), dataRaw.ts, 0.5,30,0);   
[data.imaging_frames(1,:), data.imaging_frames(2,:)] = Timing_onset_offset(dataRaw.event(14,:), dataRaw.ts, 0.5,30,0);   

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

data_signals = fieldnames(data);
for i=1:length(data_signals)
 if isempty(data.(data_signals{i}))  %remove empty fields
    data= rmfield(data, (data_signals{i}));
 end
end



%assignedResp = [1 1 1 2 2 2];
corrChoice = [-1 -1 -1 NaN 1 1 1];

% Direction - L - -1, R - 1

tasteStimEvent = []; tasteStim = []; resp = [];

% extract taste stimuli
for i = [1:3 5:7]
% for i = [1:2]
   tasteStimEvent = [tasteStimEvent, eval(sprintf('data.T_%d',i))];
   tasteStim = [tasteStim, i*ones(1,length(eval(sprintf('data.T_%d',i))))];

end

[~,ind] = sort(tasteStimEvent);
tasteStim = tasteStim(ind);

%% Re-arrange to trial structure
%it also means that trials without central sp licks are not extracted

x=2;

trialStart = data.centSp(1); trialEnd = data.centSp(1) + 9;

trial(1).sampleTime = data.centSp(1);

leftLicks = data.LeftSp((data.LeftSp < trialEnd) & (data.LeftSp > trialStart));
rightLicks = data.RightSp((data.RightSp < trialEnd) & (data.RightSp > trialStart));

leftReward = data.L_1((data.L_1 < trialEnd) & (data.L_1 > trialStart));
rightReward = data.R_1((data.R_1 < trialEnd) & (data.R_1 > trialStart));
        
trial(1).tasteStim = tasteStim(1);

if ~isempty(leftLicks) && ~isempty(rightLicks)
            
   trial(1).choice = (leftLicks(1) -  rightLicks(1))/abs(leftLicks(1) - rightLicks(1));
   trial(1).actionOnset = min(leftLicks(1),rightLicks(1));

elseif ~isempty(leftLicks) && isempty(rightLicks)

   trial(1).choice = -1;
   trial(1).actionOnset = leftLicks(1);

elseif isempty(leftLicks) && ~isempty(rightLicks)

   trial(1).choice = 1;
   trial(1).actionOnset = rightLicks(1);

end

if ~isempty(leftReward)
    trial(1).rewardOnset = leftReward;
elseif ~isempty(rightReward)
    trial(1).rewardOnset = rightReward;
end

trial(1).corrChoice = corrChoice(tasteStim(1));
                  

data.firstLick(1) = data.centSp(1);

for i =2:length(data.centSp) 
    
    if data.centSp(i)-data.centSp(i-1) <0.5
        continue
    else

        trialStart = data.centSp(i); trialEnd = data.centSp(i) + 9;
        
        data.firstLick(x) = data.centSp(i);

        trial(x).sampleTime = data.centSp(i);
        
        leftLicks = data.LeftSp((data.LeftSp < trialEnd) & (data.LeftSp > trialStart));        
        rightLicks = data.RightSp((data.RightSp < trialEnd) & (data.RightSp > trialStart));
                
        leftReward = data.L_1((data.L_1 < trialEnd) & (data.L_1 > trialStart));
        rightReward = data.R_1((data.R_1 < trialEnd) & (data.R_1 > trialStart));
        
        trial(x).tasteStim = tasteStim(x);
        
        if ~isempty(leftLicks) && ~isempty(rightLicks)
            
           trial(x).choice = (leftLicks(1) -  rightLicks(1))/abs(leftLicks(1) - rightLicks(1));
           trial(x).actionOnset = min(leftLicks(1),rightLicks(1));
           
        elseif ~isempty(leftLicks) && isempty(rightLicks)
            
           trial(x).choice = -1;
           trial(x).actionOnset = leftLicks(1);
        
        elseif isempty(leftLicks) && ~isempty(rightLicks)
            
           trial(x).choice = 1;
           trial(x).actionOnset = rightLicks(1);
           
        end
        
        if ~isempty(leftReward) && isempty(rightReward)
            trial(x).rewardOnset = leftReward;
        elseif ~isempty(rightReward) && isempty(leftReward)
            trial(x).rewardOnset = rightReward;
        elseif ~isempty(rightReward) && ~isempty(leftReward)
            trial(x).rewardOnset = 0;
        end
        
        trial(x).corrChoice = corrChoice(tasteStim(x));
                       
        x=x+1;
    end
end



%% list trials with no lateral licks 

j=1;
lateralmiss=[];
for i=1:length(trial)
    
    A = isempty(trial(i).choice);
    
    if A
        lateralmiss(j) = i;
        j=j+1;
    end
end

%% problem trials (random event recorded by national instruments results in event in both left and right water delivery
j=1;
problemTrial = [];
for i =1:length(trial)
    if trial(i).rewardOnset == 0
        problemTrial(j) = i;
        j=j+1;
        
    end
end

%% remove trial where there is no taste ID

j=1;
noTasteID = []; 
for i =1:length(trial)
    if isempty(trial(i).tasteStim)
        noTasteID(j) = i;
        j=j+1; 
    
    end
end
%% remove these trials from trial struct

a=[problemTrial  lateralmiss noTasteID];

trial(a) = [];

%% add if it was a left or right trial based on user assigned directions
% 
% for i=1:length(trial)
%     for j=1:length(tastes)
%         if convertCharsToStrings(trial(i).TasteID(:)) == tastes{j}
%             trial(i).L_R_trial = assigned_direction(j);
%         end
%     end
% end

%% add correct or incorrect trials based on if water was delivered at lateral spout

% for i =1:length(trial)
%     if isempty(trial(i).leftReward) && isempty(trial(i).rightReward)
%         trial(i).correct_choice = 0;
%     else
%         trial(i).correct_choice = 1;
%     end
% end


%% total performance
% x=0;
% for i=1:length(trial)
%     if trial(i).correct_choice == 1
%         x=x+1;
%     end
%     trial(i).total_performance = x/i;
% end
%% performance for each taste
% x=zeros(1,length(tastes));
% trialcount=zeros(1,length(tastes));
% b=1;
% for j=1:length(tastes)
%     for i=1:length(trial)
%         if trial(i).correct_choice == 1 & convertCharsToStrings(trial(i).TasteID(:)) == tastes{j}
%             x(b)=x(b)+1;   
%         end
%         
%         if convertCharsToStrings(trial(i).TasteID(:)) == tastes{j}
%             trialcount(b)=trialcount(b)+1;
%             trial(i).((append((tastes{j}),'_performance'))) = x(b)/trialcount(b);
%             
%         else
%             trial(i).((append((tastes{j}),'_performance'))) = [];
%         end
%     end
%     b=b+1;
% end

%% bias calculation are they making more errors in one direction than the other

% Lerror = 0;
% Rerror=0;
% Lcount = 0;
% Rcount = 0;
% for i =1:length(trial)
%     if trial(i).correct_choice == 0 & trial(i).L_R_trial == 1
%         Lerror = Lerror+1;
%     end
%     if trial(i).correct_choice == 0 & trial(i).L_R_trial == 2
%         Rerror = Rerror+1;
%     end
%     if trial(i).L_R_trial == 1
%         Lcount = Lcount+1;
%     end
%     if trial(i).L_R_trial == 2
%         Rcount = Rcount+1;
%     end
%     trial(i).bias = Lerror/Lcount - Rerror/Rcount;
% end
% summary.bias = Lerror/Lcount - Rerror/Rcount;
%     

%% summary performance for each taste
% T_correct = zeros(1,length(tastes));
% T_count = zeros(1,length(tastes));
% t_performance = zeros(1,length(tastes));
% for j=1:length(tastes)
%     for i =1:length(trial)
%         if convertCharsToStrings(trial(i).TasteID) == tastes{j}
%             T_count(j) = T_count(j)+1;
%             if trial(i).correct_choice == 1
%                 T_correct(j) = T_correct(j)+1;
%             end
%         end   
%     end
% end
% for j=1:length(tastes)
%     summary.(append((tastes{j}),'_performance')) = T_correct(j)/T_count(j);
% end
% t_performance = T_correct./T_count;
% summary.ind_performance = t_performance;
% 
% 
% %% summary total performance
% %all outputs of this function are also saved in the summary struct
% summary.total_performance = trial(length(trial)).total_performance;
% summary.numTrials = length(trial);
% summary.tastes = tastes;
% summary.directions = excel_directions;
% summary.trial = trial;
% summary.unit = unit;
% summary.data = data;
% summary.problemTrial = problemTrial;
% summary.lateralmiss = lateralmiss;
% summary.noTasteID = noTasteID;
%     
end