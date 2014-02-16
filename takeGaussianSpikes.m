function data = takeGaussianSpikes(data,inputs)

        nGaussComponents = 6;
        nClusters = 4;
        ppc = nGaussComponents*3;
        for clustN = 1:nClusters
            heightIX = ppc*(clustN - 1) + 1 + ([1:nGaussComponents]-1)*3;
            muIX     = ppc*(clustN - 1) + 2 + ([1:nGaussComponents]-1)*3;
            sigIX    = ppc*(clustN - 1) + 3 + ([1:nGaussComponents]-1)*3; 
            data.spikeAvg{clustN} = gaussianComposition(...
                            inputs(heightIX),...
                            inputs(muIX),...
                            inputs(sigIX),...
                            [-data.spikeWidth:1/(data.sampleRate):data.spikeWidth]);
        end

