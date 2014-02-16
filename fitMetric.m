function totalScore = fitMetric(data)
  
    LPF = data.LPF;         
    HPF = data.HPF;       
    
    % Get filtered data
    dVdT = data.dVdT;
      
    % For each cluster, convolve spikes with raster, filter it
    fullModel = dVdT - makeResidual(data);
    lpM    = filtfilt(data.LPF.d.sosMatrix,data.LPF.d.ScaleValues,fullModel);
    fModel = filtfilt(data.HPF.d.sosMatrix,data.HPF.d.ScaleValues,lpM);

    residual = dVdT - fModel;
    
    % Take power spectra of source and residual
    [PdVdT,f] = pwelch(dVdT,data.sampleRate/2,data.sampleRate/4,[],data.sampleRate);
    [Presid,f] = pwelch(residual,data.sampleRate/2,data.sampleRate/4,[],data.sampleRate);
    FsIX = dsearchn(f,HPF.freq);
    FeIX = dsearchn(f,LPF.freq);
    % Integrate power in the band of interest
    sourcePower = sum(PdVdT(FsIX:FeIX));
    residPower  = sum(Presid(FsIX:FeIX));
    
    residScore = sourcePower/residPower;
    
    totalScore = residScore;
    
    

    
    
    
        
    