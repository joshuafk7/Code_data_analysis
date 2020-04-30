function plot_licking(summary)
for i=1:length(summary.trial)
    j = ones(1,length(summary.trial(i).centSp))*i;
    h1=scatter(summary.trial(i).centSp,j,6,'filled','r'); %central licks are red
    
     hold on
end

for i=1:length(summary.trial)
    j = ones(1,length(summary.trial(i).LeftSp))*i;
    h2=scatter(summary.trial(i).LeftSp,j,6,'filled','b'); %left licks are blue
     hold on
end
for i=1:length(summary.trial)
    j = ones(1,length(summary.trial(i).RightSp))*i;
    h3=scatter(summary.trial(i).RightSp,j,6,'filled','g'); %right licks are green
     hold on
end
legend([h1,h2,h3],'central','left','right','Location','northwest')
title('Licking')
ylabel('Trial #');
xlabel('Time (sec)')
xlim([-4 6]);
ylim([0 length(summary.trial)]);
% suptitle(append(summary.mouseID))
end