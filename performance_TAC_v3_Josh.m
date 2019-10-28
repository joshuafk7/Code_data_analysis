%% extract stimuli
function performance_TAC_v3(x,tasteID)
% This function is modified by Ke to be compatible with 4 tastant
% discrimination
tasteID=tasteID+1;
cwd = pwd;
load('data.mat');
LeftTaste=[];
for i=1:length(data.leftID)
    if ~iscell(data.leftID)
        LeftTaste=data.(data.leftID).dig;
        LeftTaste(:,2)=1;
        LeftTaste(:,3)=data.(data.leftID).line;
    else
        Left=data.(data.leftID{i}).dig;
        Left(:,2)=1; % use 1 for taste associated with left spout
        Left(:,3)=data.(data.leftID{i}).line;
        LeftTaste=[LeftTaste; Left];
        clear Left
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RightTaste=[];
for i=1:length(data.rightID)
    if ~iscell(data.rightID)
        RightTaste=data.(data.rightID).dig;
        RightTaste(:,2)=2; % use 3 for taste associated with right spout
        RightTaste(:,3)=data.(data.rightID).line;
    else
        Right=data.(data.rightID{i}).dig;
        Right(:,2)=2; % use 3 for taste associated with right spout
        Right(:,3)=data.(data.rightID{i}).line;
        RightTaste=[RightTaste; Right];
        clear Right
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stimuli         =[LeftTaste;RightTaste];
%% sort stimuli
[~, Stimuli_ind]=sort(stimuli(:,1)); % sort the event (W and NaCl) as ascending
for i=1:length (Stimuli_ind)
    st(i,:)=stimuli(Stimuli_ind(i),:);
end
clear stimuli
stimuli=st;
%%
%sumData.id=filename;
sumData.mouseID = data.mouseID;
sumData.date    = data.date   ;
sumData.stimuli=stimuli;
%%

% LeftCorrect=data.Left_Correct.dig;
% RightCorrect=data.Right_Correct.dig;
decision_onset=data.Up.dig_offset;

%% 1st, find out the trial with no licks to the centeral port
trial_noSampling=[];
% for i=1:size(st,1)
%     if isempty(find(decision_onset-st(i,1)>2 & decision_onset-st(i,1)<7 ))
%         trial_noSampling=[trial_noSampling i];
%     end
% end
for i=1:size(st,1)
    if isempty(find(decision_onset-st(i,1)>0 & decision_onset-st(i,1)<5 ))
        trial_noSampling=[trial_noSampling i];
    end
end
sumData.noSampling=trial_noSampling;
%%
st(trial_noSampling,:)=[];

if size(st,1)~= length(decision_onset)
    warning('Be careful to check the data, Probably there is something wrong')
end
trial_end=data.Down.dig;
if size(st,1)~= length(trial_end)
    warning('Be careful to check the data, Probably there is something wrong')
    if size(st,1)<length(trial_end)% it means recorded more Up/Down lateral port signal
        for i=1:size(st,1)
            ind=find(trial_end>st(i,1));
            trial_end(i)=trial_end(ind(1));
        end
    else
        st(end,:)=[];   % it means recorded less Up/Down lateral port signal (likely happended for the last trial)
    end
end
trial_noResponse=[];
LeftError.trial=[];
RightError.trial=[];
for i=1:size(st,1)
    raw{i,1}=st(i,2);
    raw{i,2}=st(i,1);
    c=spike2eventRasteandPSTH_NP(data.centSp.dig,st(i,1),100,0,ceil(trial_end(i)-st(i,1))*1000); % sampling
    raw{i,3}=c.spikeraster.times;
    c=spike2eventRasteandPSTH_NP(data.Up.dig_offset,st(i,1),100,0,ceil(trial_end(i)-st(i,1))*1000); % timestamps of Decision onset
    raw{i,4}=c.spikeraster.times;
    c=spike2eventRasteandPSTH_NP(data.LeftSp.dig,st(i,1),100,0,ceil(trial_end(i)-st(i,1))*1000); % timestamps of Left licks
    raw{i,5}=c.spikeraster.times;
    c=spike2eventRasteandPSTH_NP(data.RightSp.dig,st(i,1),100,0,ceil(trial_end(i)-st(i,1))*1000); % timestamps of right licks
    raw{i,6}=c.spikeraster.times;
    if isempty(raw{i,5}) & isempty(raw{i,6})
        trial_noResponse=[trial_noResponse i];
    end
    switch raw{i,1}
        case 1 % left choice
            if isempty(raw{i,5}) & ~isempty(raw{i,6})
                LeftError.trial=[LeftError.trial i];
            end
            if ~isempty(raw{i,5}) & ~isempty(raw{i,6})
                if raw{i,6}(1) < raw{i,5}(1)
                    LeftError.trial=[LeftError.trial i];
                end
            end
        case 2 % right choice
            if isempty(raw{i,6}) & ~isempty(raw{i,5})
                RightError.trial=[RightError.trial i];
            end
            if ~isempty(raw{i,5}) & ~isempty(raw{i,6})
                if raw{i,5}(1) < raw{i,6}(1)
                    RightError.trial=[RightError.trial i];
                end
            end
    end
    raw{i,7}=st(i,3);
end

LeftError.num=length(LeftError.trial);
RightError.num=length(RightError.trial);

%% 2nd find out trial with no decision making
sumData.noResponse=trial_noResponse;
left_trial=find(st(:,2)==1);
left_trial_decision=setdiff(left_trial,trial_noResponse);
Right_trial=find(st(:,2)==2);
Right_trial_decision=setdiff(Right_trial,trial_noResponse);
%% 3.2 Let's extract the "correct response" with first lick
trial_decision=[left_trial_decision; Right_trial_decision];
trial_decision=sort(trial_decision);
trial_decision(:,2)=st(trial_decision,2);
trial_decision(:,3)=1;
% sum all error trial
if isfield(LeftError,'trial') & isfield(RightError,'trial')
    error_all=[LeftError.trial RightError.trial];
elseif isfield(LeftError,'trial')
    error_all=LeftError.trial;
elseif isfield(RightError,'trial')
    error_all=RightError.trial;
else
    error_all=[];
end
if isempty(error_all)
else
    error_all=sort(error_all);
    for i=1:length(error_all)
        index=find(trial_decision(:,1)==error_all(i));
        trial_decision(index,3)=0;
    end
end
B=cumsum(trial_decision(:,3));
for i=1:size(trial_decision,1)
    trial_decision(i,4)=B(i)/i;
end
%%

%% save data
sumData.LeftTrialDecision=left_trial_decision;
sumData.LeftError=LeftError;
sumData.RightTrialDecision=Right_trial_decision;
sumData.RightError=RightError;
sumData.MissLateral=length(trial_noResponse)/(size(st,1));
sumData.Performance=1-(LeftError.num+RightError.num)/length([left_trial_decision;Right_trial_decision]);
sumData.Performance_bias=LeftError.num/length(left_trial_decision)-RightError.num/length(Right_trial_decision);
sumData.Trial_decision=trial_decision;
sumData.raw=raw;
sumData.raw_decision=raw(trial_decision(:,1),:);
clearvars -except sumData FinalPerformance data tasteID x cwd
% save('sumData.mat','sumData')
a=num2cell(sumData.Trial_decision);
sumData.raw_decision=[a sumData.raw_decision];
sumData.raw_decision(:,5)=[];
sumData.tasteID          = x;
Performance=cell2mat(sumData.raw_decision(:,[3,10]));
for i=1:length(tasteID)
    p(i)=sum(Performance(find(Performance(:,2)==tasteID(i)),1))/length(find(Performance(:,2)==tasteID(i)));
end
sumData.performance.p = p;
sumData.performance.x = x;
cd('D:\Behavior\Discrimination\summary');

if exist(data.mouseID)==7
    cd(['D:\Behavior\Discrimination\summary\' data.mouseID]);
    load('Summary.mat');
    f = size(summaryData,2);
    summaryData(f+1) = sumData;
    save('Summary','summaryData');
else
    mkdir(data.mouseID);
    cd(['D:\Behavior\Discrimination\summary\' data.mouseID]);
    f=1;
    summaryData = sumData;
    save('Summary','summaryData');
end


%save('sumData','sumData');
%%
cd(cwd)
fig1 = figure;
subplot(1,2,2);
sumData.Trial_decision=[sumData.Trial_decision, cell2mat(sumData.raw_decision(:,end))];
for k = 1:size(sumData.Trial_decision,1)
    switch length(tasteID)
        case 2
            if sumData.Trial_decision(k,end)==tasteID(1)
                h1 = plot([0 0],[k+0.5 k+0.5],'oc', 'MarkerSize',3,'MarkerFaceColor','c');hold on;
                for i = 1:size(sumData.raw_decision{k,6},2)
                    plot([sumData.raw_decision{k,6}(i) sumData.raw_decision{k,6}(i)],[k k+1],'c');hold on;
                end
            elseif sumData.Trial_decision(k,end)==tasteID(2)
                h2 = plot([0 0],[k+0.5 k+0.5],'ob', 'MarkerSize',3,'MarkerFaceColor','b');hold on;
                for i = 1:size(sumData.raw_decision{k,6},2)
                    plot([sumData.raw_decision{k,6}(i) sumData.raw_decision{k,6}(i)],[k k+1],'b');hold on;
                end
            end
        case 4
            if sumData.Trial_decision(k,end)==tasteID(1)
                h1 = plot([0 0],[k+0.5 k+0.5],'oc', 'MarkerSize',3,'MarkerFaceColor','c');hold on;
                for i = 1:size(sumData.raw_decision{k,6},2)
                    plot([sumData.raw_decision{k,6}(i) sumData.raw_decision{k,6}(i)],[k k+1],'c');hold on;
                end
            elseif sumData.Trial_decision(k,end)==tasteID(2)
                h2 = plot([0 0],[k+0.5 k+0.5],'ob', 'MarkerSize',3,'MarkerFaceColor','b');hold on;
                for i = 1:size(sumData.raw_decision{k,6},2)
                    plot([sumData.raw_decision{k,6}(i) sumData.raw_decision{k,6}(i)],[k k+1],'b');hold on;
                end
            elseif sumData.Trial_decision(k,end)==tasteID(3)
                h3 = plot([0 0],[k+0.5 k+0.5],'or', 'MarkerSize',3,'MarkerFaceColor','r');hold on;
                for i = 1:size(sumData.raw_decision{k,6},2)
                    plot([sumData.raw_decision{k,6}(i) sumData.raw_decision{k,6}(i)],[k k+1],'r');hold on;
                end
            elseif sumData.Trial_decision(k,end)==tasteID(4)
                h4 = plot([0 0],[k+0.5 k+0.5],'ok', 'MarkerSize',3,'MarkerFaceColor','k');hold on;
                for i = 1:size(sumData.raw_decision{k,6},2)
                    plot([sumData.raw_decision{k,6}(i) sumData.raw_decision{k,6}(i)],[k k+1],'k');hold on;
                end
            end
            
    end
    if sumData.Trial_decision(k,2)==1 && sumData.Trial_decision(k,3)==1
        
        for i = 1:size(sumData.raw_decision{k,8},2)
            plot([sumData.raw_decision{k,8}(i) sumData.raw_decision{k,8}(i)],[k k+1],'r');hold on;
        end
        
    elseif sumData.Trial_decision(k,2)==2 && sumData.Trial_decision(k,3)==1
        
        for i = 1:size(sumData.raw_decision{k,9},2)
            plot([sumData.raw_decision{k,9}(i) sumData.raw_decision{k,9}(i)],[k k+1],'b');hold on;
        end
        
    elseif sumData.Trial_decision(k,2)==1 && sumData.Trial_decision(k,3)==0
        
        for i = 1:size(sumData.raw_decision{k,9},2)
            plot([sumData.raw_decision{k,9}(i) sumData.raw_decision{k,9}(i)],[k k+1],'k');hold on;
        end
        
    elseif sumData.Trial_decision(k,2)==2 && sumData.Trial_decision(k,3)==0
        
        for i = 1:size(sumData.raw_decision{k,8},2)
            plot([sumData.raw_decision{k,8}(i) sumData.raw_decision{k,8}(i)],[k k+1],'k');hold on;
        end
        
    end
    
    
    if     sumData.Trial_decision(k,2)==2     && sumData.Trial_decision(k,3)==1
        h5 =  plot([6.5 7],[k+0.5 k+0.5],'g');hold on;
    elseif sumData.Trial_decision(k,2)==1     && sumData.Trial_decision(k,3)==1
        plot([6 6.5],[k+0.5 k+0.5],'g');hold on;
    elseif sumData.Trial_decision(k,2)==1     && sumData.Trial_decision(k,3)==0
        h6 = plot([6 6.5],[k+0.5 k+0.5],'m');hold on;
    elseif sumData.Trial_decision(k,2)==2     && sumData.Trial_decision(k,3)==0
        plot([6.5 7],[k+0.5 k+0.5],'m');hold on;
    end
    
end

xlim([-1 9]);
switch length(tasteID)
    case 2
        legend([h1,h2],x)
    case 4
        legend([h1,h2,h3,h4],x);
end
xlabel('Time(s)') ;
ylabel('# Trials');
title([data.mouseID ' - ' num2str(size(summaryData,2)) ' day - Performance:' num2str(sumData.Performance)]);
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
% f=f+1;
subplot(1,2,1);
plot(4:size(sumData.Trial_decision,1),sumData.Trial_decision(4:end,4));hold on;
plot(4:size(sumData.Trial_decision,1),sumData.Trial_decision(4:end,4),'ro','MarkerFaceColor','r');hold on;
plot([4 size(sumData.Trial_decision,1)],[0.5 0.5],'--');hold on;
ylim([0.3 1]);
ff = fit((4:size(sumData.Trial_decision,1))',sumData.Trial_decision(4:end,4),'exp2');
plot(ff,(4:size(sumData.Trial_decision,1))',sumData.Trial_decision(4:end,4));
ylabel('Performance');
xlabel('Trials');
print(fig1,'performance','-dpng')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot the performance for different tastant

figure
bar(p,'FaceColor','w','EdgeColor','k','LineWidth',1.5)
% sin = 'r(1) + r(2) * sin(theta - r(3))';
% sinfxn = inline(sin, 'r', 'theta');
% r = nlinfit(x,y,sinfxn, [1 1 0]);

ylim([0 1])
set(gca,'xticklabel',x)
set(gca,'TickDir','out')
set(gca,'FontSize',12)
box off
print('performance2','-dpng')


cd('D:\Behavior\Discrimination');