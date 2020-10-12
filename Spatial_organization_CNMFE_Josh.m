%%Spatial mapping function
function [correlationmatrix, center,sortedD] = Spatial_organization_CNMFE_Josh(neuron)
%% reshape 2d matrix into 3d matrix (x by y by neuron)
%extract center of neuron based on largest value in A and make a scatter
%plot
A = neuron.reshape(neuron.A, 2);
% B = full(neuron.A);
% C = reshape(B ,[61 109 141]);
% imagesc(C(1,:,:))
for i=1:length(A(1,1,:))
     ai = A(:,:,i);
    [row,col] = find(ai==max(ai,[],'all'));
    center(i,1) = col;
    center(i,2) = row;
    center(i,3) = i;
    ai=[];
end
Coor = neuron.show_contours(0.8); 
hold on
scatter(center(:,1), center(:,2),'filled');
 text(center(:,1)+.5, center(:,2)+.5, num2str(center(:,3)),'fontsize',9);
 title('Spatial map of neurons');
 hold off
 %% find distances between pairs of cells
 B = pdist(center);
 distmatrix = squareform(B);
 
 %% find correlations between neuron pairs using C variable
%  correlationmatrix = zeros(length(neuron.S(:,1)),length(neuron.S(:,1)));
%  for i =1:length(neuron.S(:,1))
%      for j =1:length(neuron.S(:,1))
%         spikes1=find(neuron.S(i,:));
%         a=length(spikes1);
%         spikes2=find(neuron.S(j,:));
%         overlap = intersect(spikes1,spikes2);
%         c=length(overlap);
%         correlationmatrix(i,j) = c/a;
%      end
%  end
% for i = 1:length(neuron.C(:,1))
%      temp = [];
%      z = mod(length(neuron.C(i,:)),4);
%      temp = neuron.C(i,:);
%      temp(end-z+1:end) = [];
%      temp = full(reshape(temp,[],4));
%      temp1(:,i) = sum(temp,2);
%      
%  end
% %  temp2 = temp1(:,a);
%  correlationmatrix = corr(temp1);
%  temp1=[];
 temp1 = neuron.C';
 
%  temp2 = temp1(:,a);

 correlationmatrix = corr(temp1);
 
 %% make a list with distances and correlations
 %sortedD contains distance, correlation and neurons sorted by correlation
 %parameter
 
 x=1;
 D = reshape(distmatrix, [length(distmatrix)^2,1]);
 E = reshape(correlationmatrix, [length(correlationmatrix)^2,1]);
 D(:,2) = E;
  for i =1:length(neuron.C(:,1))
     for j =1:length(neuron.C(:,1))
         
         D(x,3) = i;
         D(x,4) = j;
         x=x+1;
     end
  end
F=D(:,2);
[G,idx] = sort(F,'descend');
sortedD = D(idx,:);

end