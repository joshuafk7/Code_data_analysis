
%% load data

file = dir('*.rhd');
[~,~,~,trial,summary] = process_intan_v2_behavior_only(file.name);
cd D:\Behavior\Discrimination\summary
if exist(summary.mouseID) == 7
    cd(summary.mouseID)
else
    mkdir(summary.mouseID);
    cd(summary.mouseID);
end
if exist('total_performance.mat') == 2
    load('total_performance.mat')
    total_performance = [total_performance,summary];
else
    total_performance = summary;
end


save('total_performance.mat','total_performance')

