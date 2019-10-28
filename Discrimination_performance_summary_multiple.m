
%% summarize behavioral data from 2afc rig in 2p room
%data is saved into total performance struct
%each "stage" is a different part of the behavioral training paradigm, all
%days of the same type are aggregated into the stages
function [trial,summary, total_performance] = Discrimination_performance_summary_multiple(filename,excel_tastes,excel_directions)
% file = dir('*.rhd');
A=cd;
date = str2double(A(34:end));

cd ..

if exist('summary','dir') ==7
    cd summary
else
    mkdir('summary');
    cd summary;
end
% a = length(dir())-2;
% i=length(dir())-1
i=1;
x=0;
if exist('total_performance.mat') == 2
    load('total_performance.mat','total_performance')
    k = length(fieldnames(total_performance));
    j=1;
    for i =1:k
        for q=1:length(total_performance.(append('stage_',(num2str(i)))))
        existingdates(j) = str2double(total_performance.(append('stage_',(num2str(i))))(q).date);
        j=j+1;
        end
    end
if ismember(date,existingdates)
    trial = [];
    summary=[];
   
    return;
else
    cd ..
    cd(num2str(date))
   [~,~,~,trial,summary] = process_intan_v3_behavior_only(filename,excel_tastes,excel_directions); 
end
    

    for i =1:k
        if isequal(fieldnames(total_performance.(append('stage_',(num2str(i))))), fieldnames(summary))
            total_performance.(append('stage_',(num2str(i)))) = [total_performance.(append('stage_',(num2str(i)))),summary];
            x=1;
            break;
        
        end
    end
    if x==0
        total_performance.(append('stage_',(num2str(k+1)))) = summary;
    end
else
    cd ..
    cd(num2str(date))
    [~,~,~,trial,summary] = process_intan_v3_behavior_only(filename,excel_tastes,excel_directions);
    total_performance.(append('stage_',(num2str(1)))) = summary;
end
cd ..
cd summary
 save('total_performance.mat','total_performance')

 end
