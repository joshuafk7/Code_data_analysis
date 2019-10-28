%% behaviral analysis for 2p rig
%% Step 1: load data from Intan
% add directory
clear
addpath('C:\Users\joshu\MATLAB\Data_analysis\')
%% load data
file = dir('*.rhd');
% filename = 'D:\Behavior\Discrimination\JK057\190808\behavior\JK057_190808_162309.rhd';
dataRaw = read_Intan(file.name);
%% extract the events
thr = .5;
[data.centSp.dig,~] = Timing_onset_offset(dataRaw.analog(1,:), dataRaw.ts, thr,30,0); % get the central licks
[data.LeftSp.dig,~] = Timing_onset_offset(dataRaw.analog(2,:), dataRaw.ts, thr,30,0);
[data.RightSp.dig,~] = Timing_onset_offset(dataRaw.analog(3,:), dataRaw.ts, thr,30,0);
[data.Forward.dig,data.Forward.dig_offset]              = Timing_onset_offset(dataRaw.event(10,:), dataRaw.ts, 0.5,30,0);
[data.Backward.dig,data.Backward.dig_offset]              = Timing_onset_offset(dataRaw.event(11,:), dataRaw.ts, 0.5,30,0);
[data.Up.dig,data.Up.dig_offset]                        = Timing_onset_offset(dataRaw.event(12,:), dataRaw.ts, 0.5,30,0);
[data.Down.dig,data.Down.dig_offset]                    = Timing_onset_offset(dataRaw.event(13,:), dataRaw.ts, 0.5,30,0);
%%
%% events
A = cd;
data.mouseID  = A(28:32);
data.date     = A(34:end);
clear A;
%%
[R_1, R_1_2]            = Timing_onset_offset(dataRaw.event(8,:), dataRaw.ts, 0.5,30,0);
[L_1, L_1_2]            = Timing_onset_offset(dataRaw.event(9,:), dataRaw.ts, 0.5,30,0);
[Taste_7,Taste_7_2]     = Timing_onset_offset(dataRaw.event(7,:), dataRaw.ts, 0.5,30,0);
[Taste_6,Taste_6_2]     = Timing_onset_offset(dataRaw.event(6,:), dataRaw.ts, 0.5,30,0);
[Taste_5,Taste_5_2]     = Timing_onset_offset(dataRaw.event(5,:), dataRaw.ts, 0.5,30,0);
[Taste_4,Taste_4_2]     = Timing_onset_offset(dataRaw.event(4,:), dataRaw.ts, 0.5,30,0);
[Taste_3,Taste_3_2]     = Timing_onset_offset(dataRaw.event(3,:), dataRaw.ts, 0.5,30,0);
[Taste_2,Taste_2_2]     = Timing_onset_offset(dataRaw.event(2,:), dataRaw.ts, 0.5,30,0);
[Taste_1,Taste_1_2]     = Timing_onset_offset(dataRaw.event(1,:), dataRaw.ts, 0.5,30,0);

Taste_1 = sort([Taste_1,Taste_1_2]);
Taste_2 = sort([Taste_2,Taste_2_2]);
Taste_3 = sort([Taste_3,Taste_3_2]);
Taste_4 = sort([Taste_4,Taste_4_2]);
Taste_5 = sort([Taste_5,Taste_5_2]);
Taste_6 = sort([Taste_6,Taste_6_2]);
Taste_7 = sort([Taste_7,Taste_7_2]);

L_1     = sort([L_1,L_1_2]);
R_1     = sort([R_1,R_1_2]);

%% reorganize the data
if isempty(Taste_1)
    Taste_1 =[];
else
    Taste_1(2,:) = 0;
end

if isempty(Taste_2)
    Taste_2 =[];
else
    Taste_2(2,:) = 1;
end

if isempty(Taste_3)
    Taste_3 =[];
else
    Taste_3(2,:) = 2;
end

if isempty(Taste_4)
    Taste_4 =[];
else
    Taste_4(2,:) = 3;
end

if isempty(Taste_5)
    Taste_5 =[];
else
    Taste_5(2,:) = 4;
end
if isempty(Taste_6)
    Taste_6 =[];
else
    Taste_6(2,:) = 5;
end
if isempty(Taste_7)
    Taste_7 =[];
else
    Taste_7(2,:) = 6;
end

if isempty(L_1)
    L_1 =[];
else
    L_1(2,:) = 7;
end

if isempty(R_1)
    R_1 =[];
else
    R_1(2,:) = 8;
end
%%
Taste= [Taste_1, Taste_2, Taste_3, Taste_4, Taste_5,Taste_6,Taste_7,L_1,R_1];
[timestamps1,Idx] = sort(Taste(1,:));
%timestamps1 = sort(Taste,1);
% assignments = [0 7;1 7;2 7;4 8;5 8;6 8]; %first number is line (0-6) and second is right or left 7 and 8)
%assignments = [0 8;1 8;2 8;4 7;5 7;6 7];


for i = 1:length(Idx)
    data1(1,i) = Taste(2,Idx(i));

end
%% Josh added section

% for i = 1:length(Idx)
%     data1(1,i) = Taste(2,Idx(i));
%     data1(2,i) = Taste(2,Idx(i));
% end
% data2 = data1(:,1:2:length(data1)-20);
% % performance = zeros(7,3);
% for j=1:length(assignments(:,1))
%     performance(j,1) = 0;
%     performance(j,2) = length(find(data2(2,:)==assignments(j,1))); %total trials per line
%     for z=2:1:length(data2)
%         if data2(2,z-1) == assignments(j,1)  && data2(2, z) == assignments(j,2)
%          performance(j,1) = performance(j,1)+1; %if trial is immediately followed by a correct r or l water delivery, increment this
%                
%         end
%         
%     end
%     performance(j,3) = performance(j,1)/performance(j,2);
% end
% x=0;
% % for z=2:1:length(data2)
% %         if data2(2,z-1) == 4  && data2(2, z) == 8
% %         x = x+1;
% %                
% %         end
% % end
% z=0;
% i=0;
% for z=2:1:length(data2)
%         if data2(2,z-1) == 6  && data2(2, z) == 7
%          i=i+1; %if trial is immediately followed by a correct r or l water delivery, increment this
%                
%         end
%         
%     end
  timestamps1 = timestamps1';
%%
%%
%data1(1:2)                 =[]; % remove the first two values cause they are useless.
%timestamps1(1:2)           =[]; % remove the first two values cause they are useless.
Left_Correct               = timestamps1(data1==7);
data.Left_Correct.dig      = Left_Correct(1:2:length(Left_Correct));
Right_Correct               = timestamps1(data1==8);
data.Right_Correct.dig      = Right_Correct(1:2:length(Right_Correct));
a=unique(data1);                % get the id of stimuli
tasteID=(find(a<7));
% message1=['You train the animal with ',num2str(length(tasteID)),' tastant, Please specify the tastant for each line.'];
% message2=[num2str((tasteID(:)+1)')];
% uiwait(msgbox({message1, message2}));
x=input('Specify the tastant\n');
y=input('Specify the left-right of each tastant,1 is left and 2 is right\n');
% x = {'100 suc' '80/20','60/40','50/50','40/60','20/100','100 NaCl'};
% y = [1 2];

switch length(tasteID)
case 2
    tasteidd = {'tastant_1','tastant_2'};
    data.leftID =[];
    data.rightID =[];
    for i = 1:length(tasteidd)
        tastant = timestamps1(data1==a(i));
        data.(tasteidd{i}).dig                                       =tastant(1:2:length(tastant));
        [data.(tasteidd{i}).cent_psth_raster] = spike2eventRasteandPSTH_NP (data.centSp.dig, data.(tasteidd{i}).dig, 100, -5000, 5000);
        [data.(tasteidd{i}).right_psth_raster] = spike2eventRasteandPSTH_NP (data.RightSp.dig, data.(tasteidd{i}).dig, 100, -5000, 5000);
        [data.(tasteidd{i}).left_psth_raster] = spike2eventRasteandPSTH_NP (data.LeftSp.dig, data.(tasteidd{i}).dig, 100, -5000, 5000);
        data.(tasteidd{i}).id=x{i};
        data.(tasteidd{i}).line=tasteID(i)+1;
        data.(tasteidd{i}).LR=y(i);
            data.leftID=tasteidd(find(y==1));       
            data.rightID=tasteidd(find(y==2));

    end       
case 6
    tasteidd = {'tastant_1','tastant_2','tastant_3','tastant_4','tastant_5','tastant_6'};
    data.leftID =[];
    data.rightID =[];
    for i = 1:length(tasteidd)
        tastant = timestamps1(data1==a(i));
        data.(tasteidd{i}).dig                                       =tastant(1:2:length(tastant));
        [data.(tasteidd{i}).cent_psth_raster] = spike2eventRasteandPSTH_NP (data.centSp.dig, data.(tasteidd{i}).dig, 100, -5000, 5000);
        [data.(tasteidd{i}).right_psth_raster] = spike2eventRasteandPSTH_NP (data.RightSp.dig, data.(tasteidd{i}).dig, 100, -5000, 5000);
        [data.(tasteidd{i}).left_psth_raster] = spike2eventRasteandPSTH_NP (data.LeftSp.dig, data.(tasteidd{i}).dig, 100, -5000, 5000);
        data.(tasteidd{i}).id=x{i};
        data.(tasteidd{i}).line=tasteID(i)+1;
        data.(tasteidd{i}).LR=y(i);
            data.leftID=tasteidd(find(y==1));       
            data.rightID=tasteidd(find(y==2));

    end
    case 7
    tasteidd = {'tastant_1','tastant_2','tastant_3','tastant_4','tastant_5','tastant_6','tastant_7'};
    data.leftID =[];
    data.rightID =[];
    for i = 1:length(tasteidd)
        tastant = timestamps1(data1==a(i));
        data.(tasteidd{i}).dig                                       =tastant(1:2:length(tastant));
        [data.(tasteidd{i}).cent_psth_raster] = spike2eventRasteandPSTH_NP (data.centSp.dig, data.(tasteidd{i}).dig, 100, -5000, 5000);
        [data.(tasteidd{i}).right_psth_raster] = spike2eventRasteandPSTH_NP (data.RightSp.dig, data.(tasteidd{i}).dig, 100, -5000, 5000);
        [data.(tasteidd{i}).left_psth_raster] = spike2eventRasteandPSTH_NP (data.LeftSp.dig, data.(tasteidd{i}).dig, 100, -5000, 5000);
        data.(tasteidd{i}).id=x{i};
        data.(tasteidd{i}).line=tasteID(i)+1;
        data.(tasteidd{i}).LR=y(i);
            data.leftID=tasteidd(find(y==1));       
            data.rightID=tasteidd(find(y==2));

    end

    otherwise
        error('You have put more than 4 tastant. You need to modify the code to process it')
end

save('data.mat','data')        
performance_TAC_v3_Josh(x,tasteID)   