function evaluateSort(data)

    % Recalculate spike waveforms
    data = makeSpikeAvg(data);

    colorList = ['r','g','b','m','c','y','k'];
    clusterList = unique(data.spikeClusters);

    figure;
    for clustN = clusterList 
        waveform = data.spikeAvg{clustN};
        plot(waveform,'Color',colorList(clustN)); hold on;
    end
    axis tight;
    
    % Generate scores from mean waveform correlations
    for clustN = clusterList

            waveform = data.spikeAvg{clustN};

            crossCorr = xcorr(data.dVdT,waveform,'none');
            waveAutoCorr = xcorr(waveform,waveform,'none');
            normFactor = waveAutoCorr(length(waveform));
            stSample = length(data.dVdT) - floor(length(waveform)/2);
            enSample = stSample + length(data.dVdT) - 1;
            corrScores(clustN,:) = log10(1./abs(1 - crossCorr(stSample:enSample)./normFactor));
    
    end
    
    spikeCorrScores = corrScores(:,data.spikeSamples);
    d0 = data.fV;
    d1 = data.dVdT;
    d2 = [diff(d1);0].*data.sampleRate; 
    d3 = [diff(d2);0].*data.sampleRate; 
    %totalScores = [spikeCorrScores;dN(data.spikeSamples)';d0(data.spikeSamples)';d1(data.spikeSamples)';d2(data.spikeSamples)'];
    totalScores = spikeCorrScores;
    [coeff,scores] = pca(totalScores');
    newClusters = kmeans(scores,length(clusterList));

    figure;
    wilsonPlot(scores,data.spikeClusters);

    
    figure;
    for clustN = clusterList
        subplot(2,2,clustN);
        ixs = find(data.spikeClusters == clustN);
        totalSpike = zeros(length(waveform),1);
        spikeList = [];
        nSpikes = 0;
        for ixN = 1:length(ixs)
            ix = ixs(ixN);
            centerSample = data.spikeSamples(ix);
            samples = centerSample + [-floor(length(waveform)/2) floor(length(waveform)/2)];
            if ((min(samples) > 0) && (max(samples) < length(data.fV)))
                plot(data.dVdT(samples(1):samples(2)),'Color',colorList(clustN));
                hold on;
                spikeList(end+1,:) = data.dVdT(samples(1):samples(2));
                totalSpike = totalSpike + data.dVdT(samples(1):samples(2));
                nSpikes = nSpikes + 1;
            end
            
        end
        plot(mean(spikeList,1),'Color','k');
        plot(std(spikeList,0,1)+2*min(mean(spikeList,1)),'Color','c');
        xlim([1 size(spikeList,2)]);
        plot(xlim(),[1 1].*2*min(mean(spikeList,1)),'c--');
        %ylim([-3.5 2]);
    end
    
    figure;
    plot(data.dVdT,'k'); hold on;
    for clustN = clusterList
        ix = find(data.spikeClusters == clustN);
        scatter(data.spikeSamples(ix),data.dVdT(data.spikeSamples(ix)),'MarkerEdgeColor',colorList(clustN));
    end
    residual = makeResidual(data);
    
    plot(max(data.dVdT)*1.1 + residual,'b');