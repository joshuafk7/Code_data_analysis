%% plot total performance for each session
%different phases of training in different colors
function plot_discrimination(total_performance,labels,d)

performance_plot = figure;
b=1;
a=[];
j=2;

for j = 1:length(fieldnames(total_performance))
    a = [b:b+length(total_performance.(append('stage_',num2str(j))))-1];
    b=b+length(total_performance.(append('stage_',num2str(j))));
%     subplot(2,1,1)
    scatter(a,vertcat(total_performance.(append('stage_',num2str(j))).total_performance),'filled')
    ylim([0 1]);
    hold on
end
% yyaxis right

line([0 b], [0.5 0.5],'LineStyle','--')
xlabel('Session #')
ylabel('Performance')
title({total_performance.stage_1(1).mouseID,'Performance'})
labels(length(labels)+1) = 'Chance';
lgd=legend(labels,'Location','southwest','FontSize',8);
title(lgd,'Suc concentrations used (mM)')
%% plot bias
bias_plot = figure;
b=1;
a=[];
bias_all=[];
for j = 1:length(fieldnames(total_performance))
    a = [b:b+length(total_performance.(append('stage_',num2str(j))))-1];
    b=b+length(total_performance.(append('stage_',num2str(j))));
    bias_all =  vertcat(bias_all,total_performance.(append('stage_',num2str(j))).bias);
end
scatter(1:b-1,bias_all)
ylim([-1 1])
title({total_performance.stage_1(1).mouseID,'Bias'})
xlabel('Session #')
ylabel('Left Bias         Right Bias')


nobias = line([0 b], [0 0],'LineStyle','-');
legend(nobias,'No Bias','Location','southwest')

%% plot individual tastes
% each stage is a separate figure
b=1;
a=[];
j=1;
p=[];
ii=1;
% colors = ['b','r', 'c','m','g','k','y'];
for j = 1:length(fieldnames(total_performance))
        iii = figure;
        sgtitle(append('Stage ',num2str(j)));
        p=[];
    for i=1:length(total_performance.(append('stage_',num2str(j))))
        
        for ii =1:length(d(j).excel_tastes)
            p(i,ii) = total_performance.(append('stage_',num2str(j)))(i).(append(d(j).excel_tastes(ii),'_performance'));
        end
        c=categorical(d(j).excel_tastes,d(j).excel_tastes);
        subplot(1,length(total_performance.(append('stage_',num2str(j)))),i) %change this
        
        bar(c,p(i,:)); %,colors(j));
        set(gca,'TickLabelInterpreter','none')
       ylim([0 1]);
      
        title(append('Day ',num2str(i)))
    end
    saveas(iii,append('Individual_taste_perf_Stage_',num2str(j),'.png'))
end
%% save plots
saveas(performance_plot,'performance_plot.png')
saveas(bias_plot,'bias_plot.png')

end