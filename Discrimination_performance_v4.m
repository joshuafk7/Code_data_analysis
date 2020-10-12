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
else 
    dates = [];
    total_perf = [];
end


cd ..
for j =1:length(b.Date)
    c = num2str(b.Date(j));
    if exist(c) == 7 && ~ismember(str2num(c),dates)
    cd(c);
    k=1;
    p=1;
    pattern = ["S","N","M","Q"]; %tastes in excel file must start with S or N
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
%     [total_perf]=Discrimination_performance_summary_multiple_v2(file.name,excel_tastes, excel_directions,dates,total_perf);
        [~,~, ~,~,summary]=process_intan_v4_behavior_only(file,excel_tastes, excel_directions);  
        total_perf = [total_perf summary];
%     if ~isempty(summary) && ~isempty(trial) %write info to excel
%         b.Performance(j) = summary.total_performance;
%         b.Bias(j) = summary.bias;
%         b.numTrials(j) = length(trial);
%     end
%     p=1;
%     for i=1:7
%        if ~isnan((b.(append('Dir_',num2str(i)))(j))) && ~isempty(summary)
%            b.(append('Perf_',num2str(i)))(j) = summary.(append((tastes{p}),'_performance'));
%            p=p+1;
%        end
%     end
%     if ~isempty(b.indiv_Performance)
%         b.indiv_Performance = tastes;
%         
%     end
    cd ..
    clearvars excel_directions excel_tastes
    end
    
end

cd log
% writetable(b,a(1).name) %write back to excel in log folder
cd ..
cd summary
save('total_perf.mat','total_perf')

%% extract taste names in string array

% pattern = ["S","N"]; %tastes in excel file must start with S or N
% 
% for j=1:height(b)
%     k=1;
% for i =1:7
%         if startsWith((b.(append('Line_',num2str(i)))(j)),pattern)
%             f(j).excel_tastes(k) = b.(append('Line_',num2str(i)))(j); %grab tastes from excel
%             k=k+1;
%         end
% end
% 
% end
% for i=1:length(f)
% f(i).excel_tastes = string(f(i).excel_tastes);
% end
% q=length(f);
% x=2;
% d(1).excel_tastes = f(1).excel_tastes; % this struct contains string arrays with unique tastes
% for i=2:q
%         if ~strcmp(f(i).excel_tastes(1,1), f(i-1).excel_tastes(1,1)) || length(f(i).excel_tastes) ~= length(f(i-1).excel_tastes)
%         d(x).excel_tastes = f(i).excel_tastes;
%         x=x+1;
% 
%     end

% end
% %% plot 
% labels = string(unique(b.Plot_Label,'stable'));
% plot_discrimination(total_performance,labels,d) %plotting function, also saves the updated plot in summary folder