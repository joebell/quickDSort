%%
%   corrSortFast.m
%
%   This is a spike sorting utility that takes prototype spike shapes and
%   uses them as templates for matching to the waveform. It does this by
%   computing the normalized cross-correlation of the waveform with each
%   spike shape, and then finding peaks in the cross correlation that are
%   close to 1. 
%
%   Each time it computes the cross-correlation it divides the trace into
%   chunks, and sorts only the best-matched spike in each chunk. If the
%   best matched spike lowers the standard deviation of the computed
%   residual it's kept. Otherwise, it's discarded.
%
%
function [data, lastMetric] = corrSortFast(data)

    displayOn = false;  % Allow output to screen
    chunkSize = 5;      % Spike widths
    maxZeroRounds = 10; % Consecutive rounds of sorting with zero spikes before exit.
    spikeThresh = .05;  % Spikes must lower noise power more than spikeThresh*their own power
    refractoryPeriod = .002; % Don't place same spikes within (ms)
    
    % Get the filtered voltage, and use it as the current residual
    dVdT = data.dVdT; 
    cR = dVdT;
    lastData = data;
    lastMetric = -Inf;
    data.spikeClusters = [];
    data.spikeSamples  = [];

    spikeAvg = data.spikeAvg;
    spikeWidth = length(spikeAvg{1});
    
    % Retain distinction of single and double cluster lists
    singleClusterList = 1:length(spikeAvg);
    nSingleClusters = length(singleClusterList);
    clusterList = 1:length(spikeAvg);
        
    % Compute normalization factors for the cross-correlations
    for clustN = 1:length(spikeAvg)
        waveAutoCorr = xcorr(spikeAvg{clustN},spikeAvg{clustN},'none');
        normFactor(clustN) = waveAutoCorr(length(spikeAvg{clustN}));
        spikeStd(clustN) = std(spikeAvg{clustN});
    end
    
    nZeroRounds = 0;
    totalSpikesAdded = NaN; % Seed to NaN so we cross correlate on the first run
    while (nZeroRounds <= maxZeroRounds)
        
        % For each chunk, determine a correlation score with each
        % cluster waveform
        
        % Divide the waveform into non-overlapping chunks
        % Vary the phase so the boundaries aren't in the same place each
        % cycle
        Nchunks = floor(length(dVdT)/(spikeWidth*chunkSize));
        chunkN = 1:Nchunks;
        chunkPhase = randi(chunkSize*spikeWidth);
        chunkWS = chunkPhase + (chunkN-1).*spikeWidth*chunkSize + 1;
        chunkNS = chunkWS + ceil(spikeWidth/2);
        chunkWE = chunkPhase + (chunkN-0).*spikeWidth*chunkSize;
        chunkNE = chunkWE - ceil(spikeWidth/2);
        ix = find(chunkNE > length(dVdT)); chunkNE(ix) = length(dVdT);
        ix = find(chunkWE > length(dVdT)); chunkWE(ix) = length(dVdT);
        ix = find(chunkNS > chunkNE); chunkNS(ix) = []; chunkNE(ix) = [];
        ix = find(chunkWS > chunkWE); chunkWS(ix) = []; chunkWE(ix) = [];
        
        % Only cross-correlate again if we've added new spikes
        if (totalSpikesAdded ~= 0)
            % Score the whole waveforms in each cluster
            for clustN = clusterList
                % Take the cross correlation with a spike shape
                waveform = spikeAvg{clustN};
                crossCorr = xcorr(cR,waveform,'none');           
                % Just take the central part of it, normalize it
                stSample = length(cR) - floor(length(waveform)/2);
                enSample = stSample + length(cR) - 1;
                crossCorr = crossCorr(stSample:enSample)./normFactor(clustN);
                % Store it
                allCrossCorr(clustN,:) = crossCorr; 
                % Find positive peaks in the cross correlation
                [corrPeakSamp, corrPeaks] = peakFind(crossCorr,[1,0]); % Find upper peaks
                % The score is zero everywhere there isn't a peak
                % It's maximal when the cross correlation is closest to 1
                corrScore(clustN,:) = zeros(length(cR),1);
                % Normalize scores by the spike size - give bigger spikes 
                % better scores. Also, give bigger peaks bigger scores
                corrScore(clustN,corrPeakSamp) = spikeStd(clustN)*abs(dVdT(corrPeakSamp)).*log10(1./abs(1 - corrPeaks));  
                
                % Set scores within refractoryPeriod to zero
                cIX = find(data.spikeClusters == clustN);
                existingSpikeSamples = data.spikeSamples(cIX);
                if length(existingSpikeSamples) > 0
                    closestSampleIX = dsearchn(existingSpikeSamples',corrPeakSamp);
                    distances = abs(corrPeakSamp - existingSpikeSamples(closestSampleIX)');
                    zIX = find(distances < refractoryPeriod*data.sampleRate);
                    corrScore(clustN,corrPeakSamp(zIX)) = -Inf;
                end
                % disp(clustN);
            end
        end
        
        trialData = data;
        
        % Get the best score in each chunk
        for chunkN = 1:length(chunkNS)

            for clustN = clusterList
                bestScore(clustN) = max(corrScore(clustN,[chunkNS(chunkN):chunkNE(chunkN)]));
            end
            [C, bestClust] = max(bestScore);
            [C, bestSamp]    = max(corrScore(bestClust,[chunkNS(chunkN):chunkNE(chunkN)]));
            bestSample(chunkN) = bestSamp;
            bestCluster(chunkN) = bestClust;
            
            trialData.spikeSamples(end+1)  = bestSample(chunkN) + chunkNS(chunkN) - 1;
            trialData.spikeClusters(end+1) = bestCluster(chunkN);
            if displayOn
                fprintf(num2str(bestCluster(chunkN)));
            end

        end
        if displayOn
            disp(' ');
        end
        
        % Calculate a new residual
        trialResidual = makeResidual(trialData);
        
        totalSpikesAdded = 0;
        % In each chunk, see if it's better now than it was.
        % If so, keep the spikes, otherwise toss-em
        for chunkN = 1:length(chunkWS)
            trialNoise = std(trialResidual(chunkWS(chunkN):chunkWE(chunkN)));
            cRNoise    = std(cR(chunkWS(chunkN):chunkWE(chunkN)));
            if (cRNoise - trialNoise > spikeThresh*spikeStd(clustN)/(2*chunkSize))  
            	totalSpikesAdded = totalSpikesAdded + 1;
                data.spikeSamples(end+1)  = bestSample(chunkN) + chunkNS(chunkN) - 1;
                data.spikeClusters(end+1) = bestCluster(chunkN);
                
                % Trap low ISIs
%                 cIX = find(data.spikeClusters == data.spikeClusters(end));
%                 unsortedCSamp = data.spikeSamples(cIX);
%                 cSamp = sort(unsortedCSamp,'ascend');
%                 diffSamp = diff(cSamp);
%                 min(diffSamp)
                
                if displayOn
                    fprintf('|');
                end
            elseif displayOn  
                fprintf('_');
            end
        end
        if displayOn
            disp([' ',num2str(totalSpikesAdded)]);
        end
    
        % disp(totalSpikesAdded);
        % Keep the data if we added any spikes
        if (totalSpikesAdded > 0)
           nZeroRounds = 0;
           lastData = data;
           cR = makeResidual(data);
        else
           data = lastData;
           nZeroRounds = nZeroRounds + 1;
           if nZeroRounds > maxZeroRounds
               lastMetric = fitMetric(data);
               data.sortMetric = lastMetric;
           end
        end        
           
    end
    
    
    
    
  
    
    
    