function data = makeSpikeAvg(data)

    data.spikeWidth = .010; % Sec, for a half-spike
    spikeBarrier = 1;       % Spike Widths

    dVdT = data.dVdT;

    % For each cluster make an average
    clusterList = unique(data.spikeClusters);
    spikeWidthSamp = round(data.spikeWidth*data.sampleRate);
    for clustN = clusterList
        
        % Get a spike triggered average for each cluster
        ix = find(data.spikeClusters == clustN);
        sampIxs = data.spikeSamples(ix);
        allSpikes = [];
        for sampIx = sampIxs
            % Only average well separated spikes
            spikeDiffs = sampIx - data.spikeSamples(:);
            ix = find(spikeDiffs == 0);
            spikeDiffs(ix) = [];
            if (min(abs(spikeDiffs)) > 2*spikeBarrier*spikeWidthSamp)
                stSamp = sampIx - spikeWidthSamp;
                enSamp = sampIx + spikeWidthSamp;
                if ((stSamp > 0) && (enSamp < length(dVdT)))
                    allSpikes(:,end+1) = dVdT(stSamp:enSamp);
                end
            end
        end
        % But if you don't find any, use them all
        if (size(allSpikes,2) == 0)
            for sampIx = sampIxs
                stSamp = sampIx - spikeWidthSamp;
                enSamp = sampIx + spikeWidthSamp;
                if ((stSamp > 0) && (enSamp < length(dVdT)))
                    allSpikes(:,end+1) = dVdT(stSamp:enSamp);
                end
            end
        end
        
        data.spikeAvg{clustN} = mean(allSpikes,2);
             
    end