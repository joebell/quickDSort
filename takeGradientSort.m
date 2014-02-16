function data = takeGradientSort(data, inputs)

        spikeSamples = 301;
        nClusters = 4;

        for clustN = 1:nClusters
            data.spikeAvg{clustN} = inputs((clustN-1)*spikeSamples + [1:spikeSamples]);
            %plot(data.spikeAvg{clustN},'Color',pretty(clustN)); hold on;            
        end