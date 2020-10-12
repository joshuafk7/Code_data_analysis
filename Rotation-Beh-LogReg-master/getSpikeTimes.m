function neuron = getSpikeTimes(imaging_frame_time, detectedSpikes)

for i = 1:size(detectedSpikes,1)

neuron(i).spikeTimes = imaging_frame_time(detectedSpikes(i,:) ~= 0);
neuron(i).spikeVal = detectedSpikes(i,detectedSpikes(i,:) ~= 0);

end

end

