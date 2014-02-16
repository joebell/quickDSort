%%
%   handSort.m
%
%   This function is for manually checking and fixing spike sorting. It's
%   designed to allow rapid movement through a waveform. It displays a
%   filtered trace, with different spike ID's annotated in different colors
%   on the waveform. Commands:
%
%   Navigation:
%
%       H - Pan to start
%       J - Pan left
%       L - Pan right
%       ; - Pan to end
%       I - Zoom in
%       K - Zoom out
%       O - Scan right to next spike
%       U - Scan left to next spike
%       P - Seek right to next peak in residual
%       Y - Seek left to next peak in residual
%
%       1,2 - Change threshold for seeking to residual peaks
%
%       A-G - Annotate spikes in clusters a-g. If you've scanned to a
%       residual peak, then a spike is added at that peak. If you're
%       already locked onto a spike, then that spike is converted to the
%       new cluster. If you haven't locked onto a spike or a peak,
%       cross-hairs appear and allow you to click on the trace peak to add
%       as a spike location.
%
%       R - Remove spikes. If you're locked onto a spike, that one is
%       removed and you're left in Peak Seek mode. If you're not, the last
%       spike added is removed.
%
%       W - Recalculate average spike waveforms.
%  
%       S - Save the data to disk.
%       Q - Quit sorting, return the data to the calling function.
%
%       Z - Undo all the changes from this session, return to original
%       spikes.
%
%
%                           - JSB 11/2013


function handSort(saveFilename, varargin)

if nargin > 1
    caller = varargin{1};
else
    caller = [];
end

load(saveFilename);

%% Navigation Paramters
LPF =   800;         % Hz
HPF =     1;         % Hz
zoomFactor = 1.25;
translateFactor = .25;
plotWin = .5;
peakThreshFactor = 2;       % Filters out small peaks in residual for seeking
mainWindowPosition = [ 1027         588        1916         588];
spikeAvgWindowPosition = [ 1027         230         390         274];
pointViewWindowPosition = [1421         230         390         274];
axisPosition = [.05 .05 .9 .9];

% saveFilename = 'handSortedSpikes.mat';





%%
colorList = ['r','g','b','m','c','y','k'];
spikeCategories = ['A','B','C','D','E','F','G'];
centerLineHandle = [];
spikeFix = false;
peakFix  = false;
plotPosition = plotWin/2;


% Add fields if necessary
if ~isfield(data,'spikeSamples');
    data.spikeSamples = [];
end
if ~isfield(data,'spikeClusters');
    data.spikeClusters = [];
end


if ~isfield(data,'dVdT')
    % Filter data, store
    data.LPF.freq = LPF;
    data.LPF.h = fdesign.lowpass('N,F3dB',4,data.LPF.freq/(data.sampleRate/2));
    data.LPF.d = design(data.LPF.h,'butter');
    data.HPF.freq = HPF;
    data.HPF.h = fdesign.highpass('N,F3dB',4,data.HPF.freq/(data.sampleRate/2));
    data.HPF.d = design(data.HPF.h,'butter');
    lpV = filtfilt(data.LPF.d.sosMatrix,data.LPF.d.ScaleValues,data.V);
    data.fV = filtfilt(data.HPF.d.sosMatrix,data.HPF.d.ScaleValues,lpV);
    dVdT = [diff(data.fV);0].*data.sampleRate;
    data.dVdT = dVdT;
else
    dVdT = data.dVdT;
end

% Make an average spikes waveform for calculating the residuals
if ~isfield(data,'spikeAvg')
    data = makeSpikeAvg(data);
end


% Calculate a timebase, and find all the peaks in the residual.
time = (1:length(dVdT))./data.sampleRate;
residual = []; peakTimes = []; peakHeights = []; peakSamples = [];
getResidualPeaks(data)

% Make the figures
pointFig = figure();
plotFig = figure;
avgFig = figure;
% Do the initial plots
rePlot();
plotAvgSpikes();
pointViewPlot();

% Save a copy of the data before any sorting changes
originalData = data;

    function keyPress(H,E)
        
        global selectedPointIX;
        
        % E.Key
        clustNum = 0;
        switch E.Key
            case '2'
                peakThreshFactor = peakThreshFactor*1.1;
                rePlot();
            case '1'
                peakThreshFactor = peakThreshFactor/1.1;
                rePlot();
            case 'l'
            	plotPosition = plotPosition + plotWin*translateFactor;
                xlim(plotPosition + plotWin.*[-.5 .5]);
                if spikeFix
                    delete(centerLineHandle);
                end
                if peakFix
                    delete(centerLineHandle);
                end
                spikeFix = false;
                peakFix = false;
            case 'semicolon'
            	plotPosition = time(end) - plotWin/2;
                xlim(plotPosition + plotWin.*[-.5 .5]);
                if spikeFix
                    delete(centerLineHandle);
                end
                if peakFix
                    delete(centerLineHandle);
                end
                spikeFix = false;
                peakFix = false;
            case 'h'
            	plotPosition = time(1) + plotWin/2;
                xlim(plotPosition + plotWin.*[-.5 .5]);
                if spikeFix
                    delete(centerLineHandle);
                end
                if peakFix
                    delete(centerLineHandle);
                end
                spikeFix = false;
                peakFix = false;    
            case 'j'
                plotPosition = plotPosition - plotWin*translateFactor;
                xlim(plotPosition + plotWin.*[-.5 .5]);
                if spikeFix
                    delete(centerLineHandle);
                end
                if peakFix
                    delete(centerLineHandle);
                end
                spikeFix = false;
                peakFix = false;
            case 'i'
                plotWin = plotWin/zoomFactor;
                xlim(plotPosition + plotWin.*[-.5 .5]);
            case 'k'
                plotWin = plotWin*zoomFactor;
                xlim(plotPosition + plotWin.*[-.5 .5]);
            case 'comma'
                plotWin = time(end)-time(1);
                plotPosition = plotWin/2;
                if spikeFix
                    delete(centerLineHandle);
                end
                if peakFix
                    delete(centerLineHandle);
                end
                spikeFix = false;
                peakFix = false;
                xlim(plotPosition + plotWin.*[-.5 .5]);
            % Seek right to the next spike
            case 'o'
                nearness = 1./(time(data.spikeSamples) - plotPosition);
                [C,ix] = max(nearness.*isfinite(nearness));
                if (C > 0)
                    plotPosition = time(data.spikeSamples(ix));
                    xlim(plotPosition + plotWin.*[-.5 .5]);
                    if spikeFix
                        delete(centerLineHandle);
                    end
                    if peakFix
                        delete(centerLineHandle);
                    end
                    spikeFix = true;
                    peakFix = false;
                    centerLineHandle = plot(plotPosition.*[1 1],ylim(),'k--');
                    selectedPointIX = ix;
                    % Update the point view
                    pointViewPlot();
                end                
            case 'u'
                nearness = 1./(plotPosition - time(data.spikeSamples));
                [C,ix] = max(nearness.*isfinite(nearness));
                if (C > 0)
                    plotPosition = time(data.spikeSamples(ix));
                    xlim(plotPosition + plotWin.*[-.5 .5]);
                    if spikeFix
                        delete(centerLineHandle);
                    end
                    if peakFix
                        delete(centerLineHandle);
                    end
                    spikeFix = true;
                    peakFix = false;
                    centerLineHandle = plot(plotPosition.*[1 1],ylim(),'k--');
                    selectedPointIX = ix;
                    % Update the point view
                    pointViewPlot();
                end  
            case 'p'
                nearness = 1./(peakTimes - plotPosition);
                [C,ix] = max(nearness.*isfinite(nearness));
                if (C > 0)
                    plotPosition = peakTimes(ix);
                    xlim(plotPosition + plotWin.*[-.5 .5]);
                    if spikeFix
                        delete(centerLineHandle);
                    end
                    if peakFix
                        delete(centerLineHandle);
                    end
                    spikeFix = false;
                    peakFix = true;
                    centerLineHandle = plot(plotPosition.*[1 1],ylim(),'k:');
                end                 
            case 'y'
                nearness = 1./(plotPosition - peakTimes);
                [C,ix] = max(nearness.*isfinite(nearness));
                if (C > 0)
                    plotPosition = peakTimes(ix);
                    xlim(plotPosition + plotWin.*[-.5 .5]);
                    if spikeFix
                        delete(centerLineHandle);
                    end
                    if peakFix
                        delete(centerLineHandle);
                    end
                    spikeFix = false;
                    peakFix = true;
                    centerLineHandle = plot(plotPosition.*[1 1],ylim(),'k:');
                end
            case 'a'
                clustNum = 1;
            case 'b'
                clustNum = 2;
            case 'c'
                clustNum = 3;
            case 'd'
                clustNum = 4;
            case 'e'
                clustNum = 5;
            case 'f'
                clustNum = 6;
            case 'g'
                clustNum = 7;
            % Remove a spike    
            case 'r'
                if spikeFix
                    delete(centerLineHandle);
                    spikeFix = false;
                    remIX = dsearchn(time(data.spikeSamples)',plotPosition);
                    disp(['Removed ',spikeCategories(data.spikeClusters(remIX))]);
                    remSpike(remIX);
                    peakFix = true;
                    centerLineHandle = plot(plotPosition.*[1 1],ylim(),'k:');
                elseif length(data.spikeSamples) > 0
                    disp(['Removed last ',spikeCategories(data.spikeClusters(end))]);
                    remSpike(length(data.spikeSamples)); 
                end
                rePlot();
            case 's'
                save(saveFilename,'data');
                disp(['Saved to: ',saveFilename]);
            % Update prototype display    
            case 'w'
                updateAvg();
            % Revert to original data
            case 'z'
                if strcmp(input('Are you sure you want to revert to the original sort? (y/n) ','s'),'y')
                    data = originalData;
                    plotAvgSpikes();
                    rePlot();
                end
            case 'q'
                %disp('Quit, data returned to caller.');
                close(avgFig);
                close(plotFig);
                close(pointFig);
                data.handSorted = true;
                if ~isempty(caller)
                    if caller.processFlag
                        caller.n();
                    end
                end
            otherwise
        end
        
        
        % Add or convert a spike
        if clustNum > 0
            % Convert a spike ID, stay in spikeFix mode
            if spikeFix
                spikeIX = dsearchn(time(data.spikeSamples)',plotPosition);
                disp(['Changed ',spikeCategories(data.spikeClusters(spikeIX)),...
                    ' -> ',spikeCategories(clustNum)]);
                data.spikeClusters(spikeIX) = clustNum;
                rePlot();
            % Add a spike, go to spike mode
            elseif peakFix
                peakIX = dsearchn(peakTimes,plotPosition);
                disp(['Added ',spikeCategories(clustNum)]);
                data.spikeSamples(end+1) = peakSamples(peakIX);
                data.spikeClusters(end+1) = clustNum;
                plotPosition = time(peakSamples(peakIX));
                spikeFix = true;
                peakFix = false;
                rePlot();
            else
                [x,y] = ginput(1);
                % Normalize for equal screen pixels
                aspectRatio = daspect();
                % Use peaks in raw trace, instead of residual
                [peakSamples, peakHeights] = peakFind(dVdT,[1,1]);
                ix = dsearchn([time(peakSamples)', peakHeights./aspectRatio(2)],[x,y./aspectRatio(2)]);
                addSpike(peakSamples(ix),clustNum);
                disp(['Added ',spikeCategories(clustNum)]);
                plotPosition = time(peakSamples(ix));
                spikeFix = true;
                peakFix = false;
                centerLineHandle = plot(0,0); % Make an object to delete
                rePlot();
            end
        end
    end

    function rePlot()
        
        % Update the point view
        pointViewPlot();
        
        figure(plotFig);
        clf;
        plot(time,dVdT,'k'); hold on;
        
        % Recalculate residual peaks
        getResidualPeaks(data);
        plot([time(1) time(end)],...
            (abs(min(residual)) + max(dVdT) + peakThreshFactor*std(residual)).*[1 1],...
            'k--');
        plot([time(1) time(end)],...
            (abs(min(residual)) + max(dVdT) - peakThreshFactor*std(residual)).*[1 1],...
            'k--');
        plot(time,residual + abs(min(residual)) + max(dVdT),'Color',[0 0 1]);
        
        set(plotFig,'KeyPressFcn',@keyPress);
        clusterList = unique(data.spikeClusters);
        for clustNn = 1:length(clusterList)
            clustN = clusterList(clustNn);
            ix = find(data.spikeClusters == clustN);
            scatter(time(data.spikeSamples(ix)),dVdT(data.spikeSamples(ix)),...
                'Marker','.','MarkerEdgeColor',colorList(clustN));
        end
        set(gcf,'Position',mainWindowPosition);
        set(gca,'Position',axisPosition);
        set(gca,'Color',[.8 .8 1]);        
        xlim(plotPosition + plotWin.*[-.5 .5]);  
        
        if peakFix
            centerLineHandle = plot(plotPosition.*[1 1],ylim(),'k:');
        elseif spikeFix
            centerLineHandle = plot(plotPosition.*[1 1],ylim(),'k--');
        end
    end

    function addSpike(sample, clustN)
    	data.spikeSamples(end+1) = sample;
        data.spikeClusters(end+1) = clustN;
    end

    function remSpike(ix)
        data.spikeSamples(ix) = [];
        data.spikeClusters(ix) = [];
    end

    function plotAvgSpikes()        
        figure(avgFig);
        clf;
        plotWaveforms(data);
        set(gcf,'Position',spikeAvgWindowPosition);
        figure(plotFig);
        rePlot();
    end

    function updateAvg()
        disp(['Mean waveforms updated.']);
        data = makeSpikeAvg(data);
        plotAvgSpikes();       
    end

    function getResidualPeaks(data)
        data = makeSpikeAvg(data);
        residual = makeResidual(data);
        [peakSamples, peakHeights] = peakFind(residual,[1,1]);
        % Remove low peaks to allow fast scanning.
        ix = find(abs(peakHeights) < peakThreshFactor*std(residual));
        peakSamples(ix) = [];
        peakHeights(ix) = [];
        peakTimes = time(peakSamples)';
    end

    function pointViewPlot()
        figure(pointFig);
        clf;
        plotPointView(data);
        set(gcf,'Position',pointViewWindowPosition);
        figure(plotFig);
    end

    function drawSelectedPoint()
        
        global selectedPoint;
        global selectedPointIX;
        
        if ~isempty(selectedPointIX)
            snippetBounds = round(data.spikeSamples(selectedPointIX) + [-1 1].*data.spikeWidth*data.sampleRate);
            if ((snippetBounds(1) >= 1) && (snippetBounds(2) <= length(data.dVdT)))
                snippet = data.dVdT(snippetBounds(1):snippetBounds(2));
                spikeScore(1) = sum(abs(snippet - data.spikeAvg{1}))/sum(abs(data.spikeAvg{1}));
                spikeScore(2) = sum(abs(snippet - data.spikeAvg{2}))/sum(abs(data.spikeAvg{2}));   
                selectedPoint = scatter(spikeScore(1),spikeScore(2),'ko',...
                    'SizeData',64,'LineWidth',2); hold on;
            end
        end
    end

    function plotPointView(data)
                 
        colorList = ['r','g','b','m','c','y','k','w'];
        
        % Evaluate the sort
        clusterNumbers = unique(data.spikeClusters);
        nClusters = length(clusterNumbers);
        for spikeN = 1:length(data.spikeSamples)
            snippetBounds = round(data.spikeSamples(spikeN) + [-1 1].*data.spikeWidth*data.sampleRate);
            if ((snippetBounds(1) >= 1) && (snippetBounds(2) <= length(data.dVdT)))
                snippet = data.dVdT(snippetBounds(1):snippetBounds(2));
                for clustNn = 1:nClusters
                    clustN = clusterNumbers(clustNn);
                    spikeScore(spikeN,clustN) = sum(abs(snippet - data.spikeAvg{clustN}))/sum(abs(data.spikeAvg{clustN}));
                end
            else
                for clustNn = 1:nClusters
                    clustN = clusterNumbers(clustNn);
                    spikeScore(spikeN,clustN) = 0;
                end
            end
        end
        
        for clustNn = 1:nClusters
            clustN = clusterNumbers(clustNn);
            ix = find(data.spikeClusters == clustN);
            clusterIXs{clustN} = ix;
            points = scatter(spikeScore(ix,1),spikeScore(ix,2),4,'filled','Marker','o',...
                'MarkerEdgeColor',colorList(clustN),'MarkerFaceColor',colorList(clustN),...
                'LineWidth',3);
            hold on;
            maxLim = max([xlim(),ylim()]);
            xlim([0 maxLim]); ylim([0 maxLim]);
            %line(xlim(),ylim(),'Color','k','LineStyle','--');
            set(points,'ButtonDownFcn',{@plotClick,data,spikeScore},'UserData',clustN);
        end
        
        drawSelectedPoint();
    end
        
    function plotClick(obj, event, data, spikeScore)
        
        global selectedPoint;
        global selectedPointIX;
        
        if ~isempty(selectedPoint)
            delete(selectedPoint);
            selectedPoint = [];
            selectedPointIX = [];
        end
        
        ax = get(obj, 'Parent');
        pt = get(ax, 'CurrentPoint');
        clustN = get(obj,'UserData');
        Xclick = pt(1,1);
        Yclick = pt(1,2);
        
        ix = find(data.spikeClusters == clustN);
        score1s = spikeScore(ix,1);
        score2s = spikeScore(ix,2);
        distances = (score1s - Xclick).^2 + (score2s - Yclick).^2;
        [B,distIX] = sort(distances,'ascend');
        spikeIX = ix(distIX(1));
        disp(['Cluster ',num2str(clustN),' spike #',num2str(spikeIX)]);
        
        selectedPointIX = spikeIX;
        drawSelectedPoint();
        
        % Travel to spike
        figure(plotFig);
        plotPosition = time(data.spikeSamples(spikeIX));
        xlim(plotPosition + plotWin.*[-.5 .5]);
        if spikeFix
            delete(centerLineHandle);
        end
        if peakFix
            delete(centerLineHandle);
        end
        spikeFix = true;
        peakFix = false;
        centerLineHandle = plot(plotPosition.*[1 1],ylim(),'k--');
        
    end
        
        
end
    
    
    
