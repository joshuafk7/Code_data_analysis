cd log
a = dir('*.xlsx');
b=readtable(a(1).name);
cd ..
for j =1:length(b.Date)
    c = num2str(b.Date(j));
    cd(c);
    k=1;
    p=1;
    pattern = ["S","N"];
    for i =1:7
        if startsWith((b.(append('Line_',num2str(i)))(j)),pattern)
            excel_tastes(k) = b.(append('Line_',num2str(i)))(j);
            k=k+1;
        end
        if ~isnan((b.(append('Dir_',num2str(i)))(j)))
            excel_directions(p) = (b.(append('Dir_',num2str(i)))(j));
            p=p+1;
        end
    end
    % excel_directions = cell2mat(excel_directions)
    file = dir('*.rhd');
    
    [trial,summary, total_performance]=Discrimination_performance_summary_multiple(file.name,excel_tastes, excel_directions);  
    if ~isempty(summary) && ~isempty(trial)
        b.Performance(j) = summary.total_performance;
        b.Bias(j) = summary.bias;
        b.numTrials(j) = length(trial);
    end
    cd ..
    clearvars excel_directions excel_tastes
end

cd log
writetable(b,a(1).name)