%%
%
%  clusterStructure is a cell array including instructions for
%  hierarchical clustering. Each entry in the cell array includes the ID of
%  a cluster to split (use 0 to split all events), and the IDs to split
%  into.  For example: {[0,1,3,4],[1,1,2]} First sorts all spikes into bins
%  1,3, and 4, and then takes all the 1 spikes and splits them into 1 and
%  2.
%
%
function data = quickSort(data, clusterStructure)

    plotOn = true;
    plotColors = ['r','g','b','m','y','k','c'];

    negPeaks = true;    % Negative going spikes
    peakLimit = 2;      % SD (abs. value)
    LPF =  800;         % Hz
    HPF =    1;         % Hz
   
    % Filter data, store 
    data.LPF.freq = LPF;
    data.LPF.h = fdesign.lowpass('N,F3dB',4,data.LPF.freq/(data.sampleRate/2));
    data.LPF.d = design(data.LPF.h,'butter');
    data.HPF.freq = HPF;
    data.HPF.h = fdesign.highpass('N,F3dB',4,data.HPF.freq/(data.sampleRate/2));
    data.HPF.d = design(data.HPF.h,'butter');
    
    lpV = filtfilt(data.LPF.d.sosMatrix,data.LPF.d.ScaleValues,data.V);
    filteredV = filtfilt(data.HPF.d.sosMatrix,data.HPF.d.ScaleValues,lpV);
    dVdT = [diff(filteredV);0].*data.sampleRate;
    data.fV  = filteredV;
    data.dVdT = dVdT;
    
    % Calculate derivatives on d1
    d2 = [diff(dVdT);0].*data.sampleRate;
    d3 = [diff(d2);0].*data.sampleRate;    
    % Assemble data matrix
    dM = [dVdT,d2,d3,filteredV];
    
    % Detect peaks in V
    peakList = [];
    for n=1:size(dM,2)       
        % Z-score components of dM
        component = dM(:,n);
        component = (component - mean(component))./std(component);
        dM(:,n) = component;
    end    
    
    % Remove low peaks
    if negPeaks
        [peakList,peakHeights] = peakFind(dM(:,1),[0,1]);
        ix = find(peakHeights > -peakLimit);
    else
        [peakList,peakHeights] = peakFind(dM(:,1),[1,0]);
        ix = find(peakHeights < peakLimit);
    end   
    peakList(ix) = [];
    peakHeights(ix) = [];
   
    % Do PCA on derivatives at peak points
    descriptors = dM(peakList,:);
    [coeff, score] = princomp(descriptors);
    
    % Sort into clusters hierarchically
    for clustLevel = 1:length(clusterStructure)
        
        clustInstructions = clusterStructure{clustLevel};
        getCluster      = clustInstructions(1);
        labelsToUse     = clustInstructions(2:end);
        if (getCluster == 0)
            subIXList = 1:size(descriptors,1);
        else
            subIXList = find(T == getCluster);
        end
        subDesc = descriptors(subIXList,:);
        Tsub = subCluster(subDesc,labelsToUse);
        T(subIXList) = Tsub;
    end
    
    Nclust = length(unique(T));

    
    if plotOn
        fig1 = figure;
        set(fig1,'Position',[463   781   560   420]);
        wilsonPlot(descriptors,T);
        fig2 = figure;
        set(fig2,'Position',[463   277   560   420]);
        wilsonPlot(score,T);
    end
   

    data.spikeTimes = [];
    data.spikeSamples = [];
    data.spikeClusters = [];
    for N=1:Nclust
        ix = find(T == N);       
        data.spikeTimes =    [data.spikeTimes,peakList(ix)'./data.sampleRate];
        data.spikeSamples =  [data.spikeSamples,peakList(ix)'];
        data.spikeClusters = [data.spikeClusters,N.*ones(length(ix),1)'];
    end
    data = makeSpikeAvg(data);

 if plotOn
    fig3 = figure;
    set(fig3,'Position',[1027         756         560         420]);
    plot([1:size(dM,1)]./data.sampleRate,dVdT,'k');
    hold on;
    
    for n=1:Nclust
        ix = find(data.spikeClusters == n);
        plot(data.spikeTimes(ix),dVdT(data.spikeSamples(ix)),['o',plotColors(n)]);
    end
    xlabel('Time (s)');
    ylabel('dVdT (mV/s)');
 
     pause(1);
     close(fig1);
     close(fig2);
     close(fig3);
 end
    
    
    

        
    
    
    