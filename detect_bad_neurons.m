q=zeros(1,size(neuron.C,1));
for i = 1:length(neuron.C_raw(:,1)) %find all with C_raw values less than 1 and remove
    if max(neuron.C_raw(i,:)) < 2.5
    q(i) = 1;
    end
end
length(find(q))
b=max(neuron.C_raw,[],2);
c = find(b<2.5);
d=find(cor<.9 & cor>.8);
e = intersect(c,d);

%%
i=1;
for i =1 :length(neuron.C_raw(:,1))
    AA=[];
    AA = find(full(neuron.S(i,:)) == 0);
    CC = find(full(neuron.S(i,:)) ~= 0);
    Bb(1,i) = var(neuron.C_raw(i,AA));
    Bb(2,i) = var(neuron.C_raw(i,CC));
end
Bb(3,:)= Bb(2,:)-Bb(1,:);
histogram(Bb(3,:))
Z = find(Bb(3,:)>.1);
%%
fakeneuron = find(q);
neuron_temp = neuron;
neuron_temp.A(:, fakeneuron)=[];
neuron_temp.C(fakeneuron,:) =[];
neuron_temp.C_raw(fakeneuron,:)=[];
neuron_temp.S(fakeneuron,:) =[];
neuron_temp.ids(fakeneuron) =[];
neuron_temp.tags(fakeneuron) =[];
neuron_new = neuron_temp;
%%
for i =1:length(neuron.C_raw(:,1)) %calculate correlation between C and C_raw
R=corrcoef(neuron.C(i,:), neuron.C_raw(i,:));
cor(i) = R(2,1);
end
%%
r=[];
r = find(cor<.8); %find cells with correlation greater than .9
%%
neuron2plot = 335;
plot(neuron.C_raw(neuron2plot,:))%,'color','g')
hold on
plot(neuron.C(neuron2plot,:))
%%
fakeneuron=[];
fakeneuron = r;
neuron_temp = neuron_new;
neuron_temp.A(:, fakeneuron)=[];
neuron_temp.C(fakeneuron,:) =[];
neuron_temp.C_raw(fakeneuron,:)=[];
neuron_temp.S(fakeneuron,:) =[];
neuron_temp.ids(fakeneuron) =[];
neuron_temp.tags(fakeneuron) =[];
neuron_new2 = neuron_temp;