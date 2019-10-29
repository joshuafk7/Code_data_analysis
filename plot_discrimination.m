%% plot total performance for each session
%different phases of training in different colors
function plot_discrimination(total_performance,labels)

plotted_data = figure;
b=1;
a=[];
j=2;

for j = 1:length(fieldnames(total_performance))
    a = [b:b+length(total_performance.(append('stage_',num2str(j))))-1];
    b=b+length(total_performance.(append('stage_',num2str(j))));
    subplot(2,1,1)
    scatter(a,vertcat(total_performance.(append('stage_',num2str(j))).total_performance),'filled')
    ylim([0 1]);
    hold on
end
% yyaxis right

line([0 b], [0.5 0.5],'LineStyle','--')
% xlabel('Session #')
ylabel('Performance')
title(total_performance.stage_1(1).mouseID)
labels(length(labels)+1) = 'Chance';
lgd=legend(labels,'Location','southwest','FontSize',8);
title(lgd,'Suc concentrations used (mM)')
%% plot bias
subplot(2,1,2)
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
title('Bias')
xlabel('Session #')
ylabel('Left Bias         Right Bias')


nobias = line([0 b], [0 0],'LineStyle','-');
legend(nobias,'No Bias','Location','southwest')

saveas(plotted_data,'performance_plot')
end