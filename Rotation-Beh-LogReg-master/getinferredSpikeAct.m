function binArray = getinferredSpikeAct(spikeStruct, eventTime, win, binWidth, tLim)

binStart = win(1) : binWidth : win(2) - binWidth;
binEnd = win(1) + binWidth : binWidth : win(2);

binArray = nan(length(eventTime),length(binStart));

for k = 1:length(eventTime)

    spikeValTr = spikeStruct.spikeVal((spikeStruct.spikeTimes > eventTime(k) + win(1)) & (spikeStruct.spikeTimes <= eventTime(k) + win(2)));
    spikeTimeTr = spikeStruct.spikeTimes((spikeStruct.spikeTimes > eventTime(k) + win(1)) & (spikeStruct.spikeTimes <= eventTime(k) + win(2)));
    
    for i = 1:length(binStart)
        
        binArray(k,i) = sum(spikeValTr((spikeTimeTr > eventTime(k) + binStart(i)) & (spikeTimeTr <= eventTime(k) + binEnd(i)))); 
        
    end

    if ~isempty(tLim)    
        binArray(k,(eventTime(k) + binStart) > tLim(k)) = NaN;
    end
    
end

end