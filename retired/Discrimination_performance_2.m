%% Behavioral performance for 2AFC mixture discrimination rig
%this script uses an excel file stored in the log folder to determine for
%each day which lines were used, which tastes they contained 
% then updates the total performance file if its not already there and
% updates the total performance plot

%it also writes the performance, # of trials, and bias to the excel sheet
clear
cd log
a = dir('*.xlsx');
b=readtable(a(1).name);
cd ..
for j =1:length(b.Date)
    c = num2str(b.Date(j));
    if exist(c) == 7
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
    p=1;
    for i=1:7
       if ~isnan((b.(append('Dir_',num2str(i)))(j))) && ~isempty(summary)
           b.(append('Perf_',num2str(i)))(j) = summary.(append((tastes{p}),'_performance'));
           p=p+1;
       end
    end
%     if ~isempty(b.indiv_Performance)
%         b.indiv_Performance = tastes;
%         
%     end
    cd ..
    clearvars excel_directions excel_tastes
    end
    
end

cd log
writetable(b,a(1).name) %write back to excel in log folder
cd ..
cd summary

%% extract taste names in string array

pattern = ["S","N"]; %tastes in excel file must start with S or N

for j=1:height(b)
    k=1;
for i =1:7
        if startsWith((b.(append('Line_',num2str(i)))(j)),pattern)
            f(j).excel_tastes(k) = b.(append('Line_',num2str(i)))(j); %grab tastes from excel
            k=k+1;
        end
end

end
for i=1:length(f)
f(i).excel_tastes = string(f(i).excel_tastes);
end
q=length(f);
x=2;
d(1).excel_tastes = f(1).excel_tastes; % this struct contains string arrays with unique tastes
for i=2:q
        if ~strcmp(f(i).excel_tastes(1,1), f(i-1).excel_tastes(1,1)) || length(f(i).excel_tastes) ~= length(f(i-1).excel_tastes)
        d(x).excel_tastes = f(i).excel_tastes;
        x=x+1;

    end

end
%% plot 
labels = string(unique(b.Plot_Label,'stable'));
plot_discrimination(total_performance,labels,d) %plotting function, also saves the updated plot in summary folder