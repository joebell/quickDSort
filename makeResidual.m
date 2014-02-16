function residual = makeResidual(data)
   
    % For each cluster, convolve spikes with raster
    clusterList = unique(data.spikeClusters);
    for clustNn = 1:length(clusterList)
       clustN = clusterList(clustNn); 
       % Make a raster 
       raster = zeros(length(data.dVdT),1);
       ix = find(data.spikeClusters == clustN);
       sampIxs = data.spikeSamples(ix);
       % sampIxs
       raster(sampIxs) = 1;
       % Convolve spike shape with raster
       modelComp(:,clustN) = conv(raster,data.spikeAvg{clustN},'same');      
       
    end
    if length(clusterList) < 1
        modelComp = zeros(length(data.dVdT),1);
    end    
    
    modelV = sum(modelComp,2);
    residual = data.dVdT - modelV;