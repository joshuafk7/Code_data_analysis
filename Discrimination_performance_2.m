%% Behavioral performance for 2AFC mixture discrimination rig
%this script uses an excel file stored in the log folder to determine for
%each day which lines were used, which tastes they contained 
% then updates the total performance file if its not already there and
% updates the total performance plot

%it also writes the performance, # of trials, and bias to the excel sheet

cd log
a = dir('*.xlsx');
b=readtable(a(1).name);
cd ..
for j =1:length(b.Date)
    c = num2str(b.Date(j));
    cd(c);
    k=1;
    p=1;
    pattern = ["S","N"]; %tastes in excel file must start with S or N
    for i =1:7
        if startsWith((b.(append('Line_',num2str(i)))(j)),pattern)
            excel_tastes(k) = b.(append('Line_',num2str(i)))(j); %grab tastes from excel
            k=k+1;
        end
        if ~isnan((b.(append('Dir_',num2str(i)))(j)))
            excel_directions(p) = (b.(append('Dir_',num2str(i)))(j)); %grab direction from excel
            p=p+1;
        end
    end
    file = dir('*.rhd');
    %run behavior analysis script
    [tastes,trial,summary, total_performance]=Discrimination_performance_summary_multiple(file.name,excel_tastes, excel_directions);  
    if ~isempty(summary) && ~isempty(trial) %write info to excel
        b.Performance(j) = summary.total_performance;
        b.Bias(j) = summary.bias;
        b.numTrials(j) = length(trial);
    end
    cd ..
    clearvars excel_directions excel_tastes
end

cd log
writetable(b,a(1).name) %write back to excel in log folder
cd ..
cd summary
plot_discrimination(total_performance) %plotting function, also saves the updated plot in summary folder