%% plotting performance across days

load('total_performance.mat')

temp = fieldnames(total_performance);
tastes=temp(4:9);
for j=1:length(total_performance)
    for i=1:length(tastes)
        x(j,i) = total_performance(j).(tastes{i});
    end
end
taste = 100:-20:0;

figure;
title('JK062')

for i =1:length(x(:,1))
        
        subplot(3,4,i)
        bar(x(i,:))
        title(append('Day ',  num2str(i)) )
        
   
end