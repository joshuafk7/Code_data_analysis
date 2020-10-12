interneurons=unique(overlapping_pair_filtered2(:,2));

% figure
% hold on
% 
% imagesc(CNMF_compressed)

figure
hold on
A = neuron.reshape(neuron.A, 2);
for i=1:length(interneurons)
     ai = A(:,:,interneurons(i));
    [row,col] = find(ai==max(ai,[],'all'));
    center(i,1) = col;
    center(i,2) = row;
    center(i,3) = i;
    ai=[];
end
% Coor = neuron.show_contours(0.8); 
% hold on
% scatter(center(:,1), center(:,2),'filled');
%%
figure(9932)
imshow('JK145_mcherry_sameplane.jpeg')
hold on
for i=1:length(neurons)
    plot(neuron.Coor{i,1}(1,:),neuron.Coor{i,1}(2,:),'Color','g')
    hold on
end
title('All gCamp ROIs identified by CNMFE')
set(gcf,'Renderer','painters')

figure(9935)
imshow('JK145_mcherry_sameplane.jpeg')
hold on
for i=1:length(ROItest)
    plot(ROItest{1,i}.mnCoordinates(:,1),ROItest{1,i}.mnCoordinates(:,2),'Color','r')
    hold on
end
title('All mCherry ROIs identified by thresholding')
set(gcf,'Renderer','painters')
%%
figure(8473)
subplot(3,1,1)
interneurons=unique(overlapping_pair(:,2));
imshow('JK145_mcherry_sameplane.jpeg')
hold on
for i=1:length(interneurons)
    plot(neuron.Coor{interneurons(i),1}(1,:),neuron.Coor{interneurons(i),1}(2,:))
    hold on
end
title('All CNMF ROIs Overlapping with mCherry')
text(50,620,append(num2str(length(interneurons)),'/',num2str(length(neurons)),' putative interneurons'),'Color','white','FontSize',8);

subplot(3,1,2)
interneurons=unique(overlapping_pair_filtered(:,2));
imshow('JK145_mcherry_sameplane.jpeg')
hold on
for i=1:length(interneurons)
    plot(neuron.Coor{interneurons(i),1}(1,:),neuron.Coor{interneurons(i),1}(2,:))
    hold on
end
title('Filter out multiple gCaMP on same interneuron')
text(50,620,append(num2str(length(interneurons)),'/',num2str(length(neurons)),' putative interneurons'),'Color','white','FontSize',8);

subplot(3,1,3)
interneurons=unique(overlapping_pair_filtered2(:,2));
imshow('JK145_mcherry_sameplane.jpeg')
hold on
for i=1:length(interneurons)
    plot(neuron.Coor{interneurons(i),1}(1,:),neuron.Coor{interneurons(i),1}(2,:))
    hold on
end
title('Require at least 75% overlap')
text(50,620,append(num2str(length(interneurons)),'/',num2str(length(neurons)),' putative interneurons'),'Color','white','FontSize',8);

set(gcf,'Renderer','painters')