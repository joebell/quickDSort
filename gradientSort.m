function gradientSort(fileName, nClusters, timeLimit)

    popSize = 100;
    load(fileName);
    
    % Use the initial seed to get started
    data = makeSpikeAvg(data);
    initVector = [];
    for clustN = 1:length(data.spikeAvg)
        initVector = [initVector,data.spikeAvg{clustN}'];
    end
    initRange = 100; % mV/sec
    
    data.spikeClusters = [];
    data.spikeSamples = [];
    data.spikeWidth = .010;
    heightRange = [-2000 2000];       % mV/s
    spikeSamples = round(2*data.spikeWidth*data.sampleRate)+1;
    
    LB = heightRange(1)*ones(1,nClusters*spikeSamples);
    UB = heightRange(2)*ones(1,nClusters*spikeSamples);
          
    gaOpts = gaoptimset('PopulationSize',popSize,...
        'TimeLimit',60*60*timeLimit,'Display','iter',...
        'PopInitRange',[initVector - initRange; initVector+initRange]);
    
    xout = ga(@fitnessWrapper,nClusters*spikeSamples,[],[],[],[],...
        LB, UB,[],gaOpts);
    
    save(fileName,'data','xout');
    
    function metric = fitnessWrapper(inputs)
        
        %clf;
        for clustN = 1:nClusters
            data.spikeAvg{clustN} = inputs((clustN-1)*spikeSamples + [1:spikeSamples]);
        %    plot(data.spikeAvg{clustN},'Color',pretty(clustN)); hold on;            
        end
        %pause(.1);
        [data, thisMetric] = corrSortFast(data);
        metric = -thisMetric;
    end

end
        
        
    