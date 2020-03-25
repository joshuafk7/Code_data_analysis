for i = 1:length(neuron.C_raw(:,1)) %find all with C_raw values less than 1 and remove
    q(i) = isempty(find(neuron.C_raw(i,:) > 1));
end
fakeneuron = find(q);
neuron_temp = neuron;
neuron_temp.A(:, fakeneuron)=[];
neuron_temp.C(fakeneuron,:) =[];
neuron_temp.C_raw(fakeneuron,:)=[];
neuron_temp.S(fakeneuron,:) =[];
neuron_temp.ids(fakeneuron) =[];
neuron_temp.tags(fakeneuron) =[];
neuron_new = neuron_temp;

for i =1:length(neuron_new.C_raw(:,1)) %calculate correlation between C and C_raw
R=corrcoef(neuron_new.C(i,:), neuron_new.C_raw(i,:));
cor(i) = R(2,1);
end
r=[];
r = find(cor<.8); %find cells with correlation greater than .9

neuron2plot = 375;
plot(neuron_new.C_raw(neuron2plot,:))
hold on
plot(neuron_new.C(neuron2plot,:))

fakeneuron = r;
neuron_temp = neuron_new;
neuron_temp.A(:, fakeneuron)=[];
neuron_temp.C(fakeneuron,:) =[];
neuron_temp.C_raw(fakeneuron,:)=[];
neuron_temp.S(fakeneuron,:) =[];
neuron_temp.ids(fakeneuron) =[];
neuron_temp.tags(fakeneuron) =[];
neuron_new2 = neuron_temp;