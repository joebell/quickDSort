function plotWaveforms(data)


colorList = ['r','g','b','m','c','y','k','w'];
    if isfield(data,'spikeAvg');
        for clustN = 1:length(data.spikeAvg)
            waveform = data.spikeAvg{clustN};
            plot(waveform,'Color',colorList(clustN)); hold on;
        end
        axis tight;
    end