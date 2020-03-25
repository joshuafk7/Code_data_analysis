%% summarize behavioral data from 2afc rig in 2p room
%data is saved into total performance struct
%each "stage" is a different part of the behavioral training paradigm, all
%days of the same type are aggregated into the stages
function [total_perf] = Discrimination_performance_summary_multiple_v2(filename,excel_tastes,excel_directions,dates)
% file = dir('*.rhd');
A=cd;
date = str2double(A(end-5:end));

cd ..

if exist('summary','dir') ==7
    cd summary
else
    mkdir('summary');
    cd summary;
end
if isfile('total_perf.mat')
   load('total_perf.mat','total_perf')
   dates = zeros(length(total_perf));
   dates = vertcat(total_perf.date);
   dates = str2num(dates);


% a = length(dir())-2;
% i=length(dir())-1
i=1;
x=0;

if ~ismember(dates,date)
    cd ..
    cd(num2str(date))
   [~,~, ~,~,summary] = process_intan_v4_behavior_only(filename,excel_tastes,excel_directions); 
   total_perf = [total_perf summary];
end
    

    
else
    cd ..
    cd(num2str(date))
    [~,~,~,~,summary] = process_intan_v4_behavior_only(filename,excel_tastes,excel_directions);
    total_perf = summary;
end
cd ..
cd summary
 save('total_perf.mat','total_perf')

 end
