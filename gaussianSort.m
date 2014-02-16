function gaussianSort(fileName, nClusters, timeLimit)
    
    load(fileName);
    data.spikeWidth = .010;
    nGaussComponents = 6;
    
    % Seed the initial fit
    ampRange = [-2000, 2000];
    muRange = [-.008,.008];
    sigRange = [0, .03];
    LBpart = [ampRange(1) muRange(1) sigRange(1)];
    UBpart = [ampRange(2) muRange(2) sigRange(2)];
    LB = []; UB = [];
    for nComp = 1:nGaussComponents
        LB = [LB, LBpart];
        UB = [UB, UBpart];
    end
    
    initVector = [];
    for clustN=1:length(data.spikeAvg)
        t = [-data.spikeWidth:(1/data.sampleRate):data.spikeWidth]';
        waveform = data.spikeAvg{clustN};
        % fitOpts = fitoptions('Method','NonLinearLeastSquares');
        fitObj = fit(t, waveform, 'gauss6','Lower',LB,'Upper',UB,'StartPoint',100*randn(nGaussComponents*3,1))
        coeffs = coeffvalues(fitObj);
        initVector = [initVector,coeffs];
        subplot(2,2,clustN);
        plot(t,waveform,'Color',pretty(clustN)); hold on;
        plot(t,fitObj(t),'k');
    end
    
    return;
        
    % Clear the existing clusters
    data.spikeClusters = [];
    data.spikeSamples = [];
   
    % Make LB and UB for GA fit
    LB = []; UB = [];
    oneLRange = [ampRange(1),muRange(1),sigRange(1)];
    oneURange = [ampRange(2),muRange(2),sigRange(2)];
    for clustN = 1:nClusters
        for gaussN = 1:nGaussComponents
            LB = [LB, oneLRange];
            UB = [UB, oneURange];  
        end
    end      
    
    initRange = .005*(UB - LB);
    gaOpts = gaoptimset('PopInitRange',[initVector - initRange;initVector + initRange],'PopulationSize',10,...
        'TimeLimit',60*60*timeLimit,'Display','iter');
    
    xout = ga(@fitnessWrapper,nGaussComponents*3*nClusters,[],[],[],[],...
        LB, UB,[],gaOpts);
    
    save(fileName,'data','xout');
    
    
    
    
    function metric = fitnessWrapper(inputs)
        
        ppc = nGaussComponents*3;
        %clf;
        for clustN = 1:nClusters
            heightIX = ppc*(clustN - 1) + 1 + ([1:nGaussComponents]-1)*3;
            muIX     = ppc*(clustN - 1) + 2 + ([1:nGaussComponents]-1)*3;
            sigIX    = ppc*(clustN - 1) + 3 + ([1:nGaussComponents]-1)*3; 
            data.spikeAvg{clustN} = gaussianComposition(...
                            inputs(heightIX),...
                            inputs(muIX),...
                            inputs(sigIX),...
                            [-data.spikeWidth:1/(data.sampleRate):data.spikeWidth]);
            %plot(data.spikeAvg{clustN},'Color',pretty(clustN)); hold on;   
            %axis tight;
        end
        %pause(.1);
        [data, thisMetric] = corrSortFast(data);
        metric = -thisMetric;
    end

    

end
        
        
    